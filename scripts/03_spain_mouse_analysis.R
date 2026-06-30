# =============================================================================
# Spain mouse cohort — scRNA-seq analysis and GSEA
# =============================================================================
# Analysis outline:
#   Part 1 — Data loading, integration, and annotation
#   Part 2 — Gene expression extraction from plasma cells
#   Part 3 — Pseudobulk GSEA: BIC vs PBIC
#   Part 4 — Pseudobulk GSEA: MIC vs PBIC
# =============================================================================

library(Seurat)
library(dplyr)
library(future); plan("multicore", workers = 15)
library(reshape2)
library(ggplot2); library(patchwork)
library(ggpubr)
library(rstatix)
library(clusterProfiler)
library(ComplexHeatmap)

setwd("...")

options(future.globals.maxSize = 1000*50 * 1024^2) # 50 GB

table = function(..., useNA = 'ifany') base::table(..., useNA = useNA)

color_clusters <- c("#dcd086", "#7535de", "#7be038", "#d939ce", "#66de76", "#532b9a",
                    "#dfd338", "#af5fd3", "#b6d85d", "#656fd8", "#5aa040", "#d662b1",
                    "#62dca7", "#df3e25", "#6cdbda", "#d63b74", "#37662c", "#8c3372",
                    "#b1d8a6", "#3d2664", "#ce963b", "#576194", "#e1773a", "#81a6dc",
                    "#c43f43", "#5ba5bc", "#914b25", "#c692d2", "#888836", "#352235",
                    "#bcd2cd", "#6a242c", "#5f9679", "#da7b80", "#253726", "#dbb5c6",
                    "#605131", "#456870", "#be9779", "#8e6575")


# =============================================================================
# HELPER FUNCTION: pseudobulk collapse per sample
# =============================================================================
# Aggregates (sums) raw counts per sample, with optional per-cell normalization.
collapse_by_sample <- function(
    seu,
    group.by = "sample_id",
    assay    = "RNA",
    slot     = "counts",
    adjust   = c("none","percell_mean","percell_cpm"),
    min_cells = 1,
    verbose  = TRUE
) {
  adjust <- match.arg(adjust)
  stopifnot(group.by %in% colnames(seu@meta.data))

  X <- Seurat::GetAssayData(seu, assay = assay, slot = slot)
  if (!inherits(X, "dgCMatrix")) {
    stop("Counts are not sparse (dgCMatrix). Convert or ensure slot='counts' is sparse.")
  }

  gvec <- as.character(seu@meta.data[[group.by]])
  keep <- !is.na(gvec)
  if (!all(keep)) {
    if (verbose) message("Dropping ", sum(!keep), " cells with NA in ", group.by)
    X    <- X[, keep, drop = FALSE]
    gvec <- gvec[keep]
  }

  tab         <- table(gvec)
  good_groups <- names(tab)[tab >= min_cells]
  if (length(good_groups) == 0) stop("No groups have >= min_cells.")
  if (length(good_groups) < length(tab) && verbose)
    message("Dropping ", sum(tab < min_cells), " groups with < ", min_cells, " cells.")

  idx_list <- split(seq_along(gvec), f = gvec)[good_groups]

  # sum counts across cells within each sample
  pb_sum <- do.call(cbind, lapply(idx_list, function(ix)
    Matrix::rowSums(X[, ix, drop = FALSE])
  ))
  colnames(pb_sum) <- good_groups

  n_cells      <- vapply(idx_list, length, integer(1))
  lib_sizes_sum <- Matrix::colSums(pb_sum)
  names(n_cells) <- colnames(pb_sum)

  if (adjust == "none") {
    pb_adj <- NULL
  } else if (adjust == "percell_mean") {
    pb_adj <- Matrix::t(Matrix::t(pb_sum) / n_cells)
  } else if (adjust == "percell_cpm") {
    pb_mean  <- Matrix::t(Matrix::t(pb_sum) / n_cells)
    lib_mean <- Matrix::colSums(pb_mean)
    pb_adj   <- Matrix::t(Matrix::t(pb_mean) / pmax(lib_mean, 1)) * 1e6
  }

  list(
    pseudobulk_sum = pb_sum,
    adjusted       = pb_adj,
    n_cells        = n_cells,
    lib_sizes_sum  = lib_sizes_sum,
    groups         = colnames(pb_sum),
    assay = assay, slot = slot, group.by = group.by, adjust = adjust
  )
}


# =============================================================================
# PART 1 — Data loading and processing (run once per batch, then save/load)
# =============================================================================

# --- 1.1  Load Cell Ranger outputs and build Seurat objects ------------------
# Directories are split into BIC/PBIC, MIC, and OTHER (remaining samples)
(mydirs1 <- list.files(path = "../mouse/spain_mouse/Cell_Ranger_counts/", full.names = T, pattern = "PBIC|BIC"))
(mydirs2 <- list.files(path = "../mouse/spain_mouse/Cell_Ranger_counts/", full.names = T, pattern = "MIC"))
(mydirs  <- list.files(path = "../mouse/spain_mouse/Cell_Ranger_counts/", full.names = T))
mydirs   <- setdiff(mydirs, c(mydirs1, mydirs2))

aux_seurat <- lapply(mydirs, function(x){
  print(samplename <- gsub(".*.\\/", "", x))

  aux_gex    <- suppressMessages(Read10X(paste0(x, "/sample_filtered_feature_bc_matrix")))
  aux_seurat <- CreateSeuratObject(counts = aux_gex, min.cells = 0, project = samplename)

  aux_seurat <- PercentageFeatureSet(aux_seurat, pattern = "^mt-", col.name = "percent.mt", assay = "RNA")
  aux_seurat$disease    <- gsub("_.*.", "", aux_seurat$orig.ident)
  aux_seurat$strain     <- gsub("Healthy|MGUS|MM|[0-9]|_", "", aux_seurat$orig.ident)
  aux_seurat$subject_id <- gsub("[A-Z]|_", "", aux_seurat$orig.ident)
  aux_seurat$cell_id    <- colnames(aux_seurat)
  aux_seurat$sample_id  <- aux_seurat$orig.ident

  # unique cell barcodes across samples
  aux_seurat$auxid   <- paste0(aux_seurat$cell_id, ":", aux_seurat$sample_id)
  colnames(aux_seurat) <- aux_seurat$auxid

  return(aux_seurat)
}); names(aux_seurat) <- mydirs

# --- 1.2  Merge and integrate samples ----------------------------------------
all_merge <- merge(aux_seurat[[1]], y = aux_seurat[-1], project = "202503_scMice_JACliment")
rm("aux_seurat")

all_merge <- all_merge %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%
  IntegrateLayers(method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca")

save(all_merge, file = "mouse_OTHER.RData")
all_merge[["RNA"]] <- JoinLayers(all_merge[["RNA"]])

table(all_merge$orig.ident)

# --- 1.3  Add cell-type annotations from collaborator -----------------
# Global annotation (broad cell types)
ibon_md        <- read.csv("../mouse/spain_mouse/20250328_All_annotated.csv"); dim(ibon_md)
ibon_md$auxid  <- paste0(ibon_md$Barcode, ":", ibon_md$Sample)
rownames(ibon_md) <- ibon_md$auxid
ibon_md <- subset(ibon_md, grepl(paste(unique(all_merge@meta.data$sample_id), collapse="|"), Sample))

all_merge_sub <- subset(all_merge, auxid %in% ibon_md$auxid)
all_merge_sub <- AddMetaData(all_merge_sub,
                             metadata = ibon_md[, c("Cluster","CellType")],
                             col.name = c("cluster.ibon","cell_type.ibon"))

# Detailed T-cell annotation (atlas-based)
ibon_md        <- read.csv("../mouse/spain_mouse/20250319_Tcells_cells_ATLAS_annotated.csv"); dim(ibon_md)
ibon_md$auxid  <- paste0(ibon_md$Barcode, ":", ibon_md$Sample)
rownames(ibon_md) <- ibon_md$auxid
ibon_md <- subset(ibon_md, grepl(paste(unique(all_merge@meta.data$sample_id), collapse="|"), Sample))

all_merge_sub <- AddMetaData(all_merge_sub,
                             metadata = ibon_md[, "Celltype", drop = F],
                             col.name = "cell_type_det.ibon")

# Combine broad + detailed T-cell labels into a single final annotation
all_merge_sub$cell_type_final.ibon <- ifelse(
  all_merge_sub$cell_type.ibon %in% c("Tcells CD8","Tcells CD4"),
  all_merge_sub$cell_type_det.ibon, all_merge_sub$cell_type.ibon
)
all_merge_sub$cell_type_final.ibon <- ifelse(
  is.na(all_merge_sub$cell_type_final.ibon),
  all_merge_sub$cell_type.ibon, all_merge_sub$cell_type_final.ibon
)
table(all_merge_sub$cell_type_final.ibon)

# --- 1.4  Remove Immunoglobulin genes (confound plasma-cell analysis) --------
ig_ids        <- grep("^Igf", grep("^Ig|Jchain", rownames(all_merge_sub), value = T), invert = T, value = T)
all_merge_sub <- subset(all_merge_sub, features = setdiff(rownames(all_merge_sub), ig_ids))

save(all_merge_sub, file = "all_merge_sub_OTHER.20250328.RData")


# =============================================================================
# PART 2 — Gene expression
# =============================================================================

# Load the three processed Seurat objects (BIC/PBIC, MIC, OTHER)
load(".../all_merge_sub_BIC.20250328.RData");   bic   <- all_merge_sub
load(".../all_merge_sub_MIC.20250328.RData");   mic   <- all_merge_sub
load(".../all_merge_sub_OTHER.20250328.RData"); other <- all_merge_sub

# Genes of interest (key myeloma / plasma-cell markers)
targets      <- c("Trp53","Tnfrsf17","Myc","Nsd2","Whsc1","Mmset")
rn           <- rownames(bic[["RNA"]])
genes_present <- intersect(targets, rn)
if (length(genes_present) == 0) stop("None of the target symbols found in RNA rownames.")

# Helper: extract normalized expression + metadata for each dataset
extract_expr <- function(seu, genes) {
  mat    <- GetAssayData(seu, assay = "RNA", layer = "data")
  expr_t <- t(as.data.frame(mat[genes, , drop = FALSE]))
  pdata  <- as.data.frame(seu@meta.data)
  merge(pdata, expr_t, by = "row.names")
}

df.bic   <- extract_expr(bic,   genes_present)
df.mic   <- extract_expr(mic,   genes_present)
df.other <- extract_expr(other, genes_present)

all <- rbind(df.bic, df.mic, df.other)

write.table(all,
  ".../all_bic_pbic_mic_other_mouse_22OCT25.txt",
  col.names = T, row.names = F, quote = F, sep = "\t")

# Reload and filter strains (remove BCMO and empty)
all   <- read.delim(".../all_bic_pbic_mic_other_mouse_22OCT25.txt", stringsAsFactors = F)
all.1 <- all[!all$strain %in% c("BCMO",""), ]
table(all.1$strain)

# Focus on malignant plasma cells (B cells_PC)
all.plasma <- all.1[all.1$cell_type_final.ibon %in% c("B cells_PC"), ]
table(all.plasma$orig.ident)

# Average expression per sample
expr_T1 <- all.plasma[, c(2, 16:19)]
info    <- unique(all.plasma[, c(2, 6, 7)]); colnames(info)[1] <- "sample"

NEW <- NULL
for (i in seq_along(unique(expr_T1$orig.ident))) {
  SAMPLE <- unique(expr_T1$orig.ident)[i]
  SUB    <- expr_T1[expr_T1$orig.ident == SAMPLE, ]
  new    <- colMeans(SUB[, 2:5])
  NEW    <- rbind(NEW, c(SAMPLE, new))
}
colnames(NEW)[1] <- "sample"
NEW <- as.data.frame(NEW)
for (i in 2:ncol(NEW)) NEW[, i] <- as.numeric(as.character(NEW[, i]))

new.1 <- merge(NEW, info, by = "sample")

# Boxplots and pairwise Wilcoxon tests per gene per strain
boxplot(new.1$Trp53 ~ new.1$strain)
stripchart(new.1$Trp53 ~ new.1$strain,
           method = "jitter", pch = 19, cex = 0.8, col = "dodgerblue", vertical = TRUE, add = TRUE)
pairwise_wilcox_test(new.1, Trp53 ~ strain)

boxplot(new.1$Tnfrsf17.new ~ new.1$strain)
stripchart(new.1$Tnfrsf17.new ~ new.1$strain,
           method = "jitter", pch = 19, cex = 0.8, col = "dodgerblue", vertical = TRUE, add = TRUE)
pairwise_wilcox_test(new.1, Tnfrsf17.new ~ strain)


# =============================================================================
# PART 3 — GSEA on pseudobulk DEGs: PBIC vs BIC
# =============================================================================

library(edgeR)
library(DESeq2)
library(fgsea)
library(msigdbr)
library(EnhancedVolcano)

# Load BIC dataset
load(".../all_merge_sub_BIC.20250328.RData")
bic <- all_merge_sub

# Subset to malignant plasma cells
bic.plasma <- subset(bic, subset = cell_type_final.ibon == "B cells_PC")

# Single-cell Wilcoxon DE (PBIC vs BIC)
Idents(bic.plasma) <- "cell_type_final.ibon"
de_plasma <- FindMarkers(
  bic.plasma,
  ident.1 = "PBIC", ident.2 = "BIC",
  group.by = "strain",
  test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1
)
de_plasma[rownames(de_plasma) == "Tnfrsf17", ]

# --- Pseudobulk collapse: aggregate counts per sample -----------------------
bic.plasma.coll    <- collapse_by_sample(bic.plasma, group.by = "sample_id", adjust = "percell_mean")
bic.plasma.coll.md <- as.data.frame(bic.plasma@meta.data)
bic.plasma.coll    <- append(bic.plasma.coll, list(bic.plasma.coll.md))

# Build sample-level expression and metadata tables
pb.exp     <- bic.plasma.coll$pseudobulk_sum
GENES      <- rownames(pb.exp)

pb.md      <- bic.plasma.coll[[10]]
sample.md  <- unique(pb.md[, c("sample_id","strain","disease")])

pb.exp.filt <- pb.exp[rownames(pb.exp) %in% GENES, ]
pb.deg      <- merge(sample.md, t(pb.exp.filt), by.x = "sample_id", by.y = "row.names", all.x = T)
dim(pb.deg)

# DESeq2 pseudobulk differential expression
expr   <- t(pb.deg[, 4:ncol(pb.deg)])
colnames(expr) <- pb.deg$sample_id
pdata  <- pb.deg[, 1:3]; rownames(pdata) <- pdata$sample_id; pdata <- pdata[, -1]

dds  <- DESeqDataSetFromMatrix(countData = expr, colData = pdata, design = ~ strain)
keep <- rowSums(counts(dds) >= 10) >= ceiling(0.2 * ncol(dds))
dds  <- dds[keep, ]
dds  <- DESeq(dds)
res  <- results(dds, contrast = c("strain","PBIC","BIC"))

res[rownames(res) == "Tnfrsf17", ]

EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'padj',
                col = c("grey30","forestgreen","grey30","red2"))

# --- GSEA using DESeq2 Wald stat as ranking ---------------------------------
ranks <- res$stat[!is.na(res$stat)]
ranks[ranks == Inf] <- 999
names(ranks) <- rownames(res)[!is.na(res$stat)]

# Gene sets: MSigDB Hallmarks + CIN70 + plasma-cell identity
CIN70_mouse <- c(
  "Tpx2","Prc1","Foxm1","Cdk1","Cfap20","Tgif2","Mcm2","H2afv","Top2a","Pcna",
  "Ube2c","Melk","Trip13","Ncapd2","Mcm7","Rnaseh2a","Rad51ap1","Kif20a","Cdc45",
  "Mad2l1","Espl1","Ccnb2","Fen1","Ttk","Cct5","Rfc4","Atad2","Ckap5","Nup205",
  "Cdc20","Cks2","Rrm2","Elavl1","Ccnb1","Rrm1","Aurkb","Msh6","Ezh2","Ctps1",
  "Dkc1","Oip5","Cdca8","Pttg1","Cep55","H2afx","Cmas","Ncaph","Mcm10","Lsm4",
  "Mtbp","Asf1b","Zwint","Pbk","Shcbp1","Cdca3","Ect2","Cdc6","Ung","Mtch2",
  "Rad21","Actl6a","Gpi1","Ddx39a","Srsf2","Hdgf","Nxt1","Nek2","Dhcr7","Aurka",
  "Ndufab1","Fam83d","Kif4"
)

plasma_cell <- c(
  "Abcb9","Ampd1","Ankrd28","Arsa","Arx","Bhlhe41","Brsk1","Cadm1","Cadps2","Ccr10","Chpf","Cibar2",
  "Cldn14","Cldn3","Creld2","Cst6","Cthrc1","Derl3","Dnaaf1","Dnajb9","Dpep1","Ell2","Erlec1","Fcrl5",
  "Fkbp11","Fkbp2","Gmppb","Gprc5d","Hdlbp","Herpud1","Hid1","Hm13","Hsp90b1","Igf1","Il5ra","Itga8",
  "Jchain","Jsrp1","Kcnma1","Kcnn3","Kdelr1","Kdelr2","Kdelr3","Krtcap2","Lman1","Lman2","Lmf1","Lmtk3",
  "Lypd6b","Manea","Manf","Mei1","Mixl1","Moxd1","Mydgf","Myl2","Mzb1","Nxpe4","Paip2b","Parm1","Pdia4",
  "Pdia6","Perp","Pla2g2d","Plaat2","Prdx4","Prss16","Rabac1","Rasgrp3","Scnn1b","Sdc1","Sdf2l1",
  "Sec11c","Selenos","Slamf7","Spaca3","Spag4","Spats2","Spcs2","Spcs3","Ssr3","Ssr4","Tgfbr3l","Tmed9",
  "Tnfrsf17","Tubb2b","Txndc11","Txndc15","Txndc5","Ube2j1","Wnt10a","Xbp1"
)

m_df     <- msigdbr(species = "Mus musculus", category = "H")
pathways <- split(m_df$gene_symbol, m_df$gs_name)
sets     <- c(pathways, list(CIN70 = CIN70_mouse, PLASMA_CELL = plasma_cell))

fg  <- fgsea(pathways = sets, stats = ranks, nperm = 10000)
res_bic <- fg[order(fg$padj), c("pathway","NES","padj")]

write.table(as.data.frame(res_bic),
  ".../Malignant_sc_plasma_cells_GSEA_Mouse_PBIC_vs_BIC.txt",
  col.names = T, row.names = F, quote = F, sep = "\t")

# Enrichment plot for the plasma-cell identity gene set
plotEnrichment(pathway = sets[["PLASMA_CELL"]], stats = ranks) +
  labs(title = "PLASMA CELL Identity", subtitle = "PBIC vs BIC") +
  theme_minimal()

# Bar plot of significant pathways (padj < 0.1)
sig_bic <- as.data.frame(res_bic[res_bic$padj < 0.1, ])

df_bic <- sig_bic %>%
  dplyr::mutate(
    direction = ifelse(NES > 0, "PBIC Enriched", "BIC Enriched"),
    label = gsub("^HALLMARK_", "", pathway),
    label = gsub("_", " ", label),
    label = gsub("REACTIVE OXIGEN", "REACTIVE OXYGEN", label, ignore.case = TRUE)
  ) %>%
  dplyr::arrange(desc(abs(NES))) %>%
  dplyr::mutate(label = factor(label, levels = rev(unique(label[order(NES)]))))

ggplot(df_bic, aes(x = label, y = NES, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = c("PBIC Enriched" = "#E66B6B", "BIC Enriched" = "#79C4D6")) +
  coord_flip() +
  labs(title = "Malignant PC PBIC vs BIC", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5))


# =============================================================================
# PART 4 — GSEA on pseudobulk DEGs: PBIC vs MIC
# =============================================================================

load(".../all_merge_sub_MIC.20250328.RData")
mic <- all_merge_sub

# Subset plasma cells from each dataset
mic.plasma  <- subset(mic,       subset = cell_type_final.ibon == "B cells_PC")
pbic.plasma <- subset(bic.plasma, subset = strain == "PBIC")

table(mic.plasma@meta.data$sample_id)

# Merge MIC and PBIC plasma cells into a single object for cross-dataset comparison
DefaultAssay(mic.plasma)  <- "RNA"
DefaultAssay(pbic.plasma) <- "RNA"
mic.plasma$dataset  <- "mic"
pbic.plasma$dataset <- "pbic"

merged <- merge(mic.plasma, y = pbic.plasma,
                add.cell.ids = c("mic.plasma","pbic.plasma"),
                project = "merged", merge.data = FALSE)
merged <- JoinLayers(merged, assay = "RNA")

# Remove a problematic sample (MM_MIC7)
merged1 <- subset(merged, subset = sample_id != "MM_MIC7")

merged1 <- NormalizeData(merged1)
merged1 <- ScaleData(merged1, features = rownames(merged1))
merged1 <- RunPCA(merged1, features = VariableFeatures(merged1))

# Single-cell Wilcoxon DE (PBIC vs MIC)
Idents(merged1) <- "cell_type_final.ibon"
de_plasma_mic <- FindMarkers(
  merged1,
  ident.1 = "PBIC", ident.2 = "MIC",
  group.by = "strain",
  test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1
)
de_plasma_mic[rownames(de_plasma_mic) == "Tnfrsf17", ]

# --- Pseudobulk collapse and DESeq2 -----------------------------------------
merged.coll    <- collapse_by_sample(merged, group.by = "sample_id", adjust = "percell_mean")
merged.coll.md <- as.data.frame(merged@meta.data)
merged.coll    <- append(merged.coll, list(merged.coll.md))

pb.exp.mic    <- merged.coll$pseudobulk_sum
pb.md.mic     <- merged.coll[[10]]
sample.md.mic <- unique(pb.md.mic[, c("sample_id","strain","disease")])

pb.deg.mic <- merge(sample.md.mic, t(pb.exp.mic), by.x = "sample_id", by.y = "row.names", all.x = T)

expr.mic <- t(pb.deg.mic[, 4:ncol(pb.deg.mic)])
colnames(expr.mic) <- pb.deg.mic$sample_id

# Replace NAs with 0 (absent genes in some samples)
expr.mic[is.na(expr.mic)] <- 0

pdata.mic <- pb.deg.mic[, 1:3]; rownames(pdata.mic) <- pdata.mic$sample_id; pdata.mic <- pdata.mic[, -1]

dds.mic <- DESeqDataSetFromMatrix(countData = expr.mic, colData = pdata.mic, design = ~ strain)
dds.mic <- dds.mic[rowSums(counts(dds.mic)) > 0, ]
dds.mic <- DESeq(dds.mic, sfType = "poscounts")
keep    <- rowSums(counts(dds.mic) >= 10) >= ceiling(0.2 * ncol(dds.mic))
dds.mic <- dds.mic[keep, ]
dds.mic <- DESeq(dds.mic)
res.mic <- results(dds.mic, contrast = c("strain","PBIC","MIC"))

res.mic[rownames(res.mic) == "Tnfrsf17", ]

# EnhancedVolcano(res.mic, lab = rownames(res.mic), x = 'log2FoldChange', y = 'padj',
#                 col = c("grey30","forestgreen","grey30","red2"))


ranks.mic <- de_plasma_mic$stat[!is.na(de_plasma_mic$stat)]
ranks.mic[ranks.mic == Inf] <- 999
names(ranks.mic) <- rownames(de_plasma_mic)[!is.na(de_plasma_mic$stat)]

fg.mic  <- fgsea(pathways = sets, stats = ranks.mic, nperm = 10000)
res_mic <- fg.mic[order(fg.mic$padj), c("pathway","NES","padj")]

write.table(as.data.frame(res_mic),
  ".../Malignant_sc_plasma_cells_GSEA_Mouse_PBIC_vs_MIC.txt",
  col.names = T, row.names = F, quote = F, sep = "\t")

plotEnrichment(pathway = sets[["PLASMA_CELL"]], stats = ranks.mic) +
  labs(title = "PLASMA CELL Identity", subtitle = "PBIC vs MIC") +
  theme_minimal()

# Bar plot of significant pathways (padj < 0.11)
sig_mic <- as.data.frame(res_mic[res_mic$padj < 0.11, ])

df_mic <- sig_mic %>%
  dplyr::mutate(
    direction = ifelse(NES > 0, "PBIC Enriched", "MIC Enriched"),
    label = gsub("^HALLMARK_", "", pathway),
    label = gsub("_", " ", label),
    label = gsub("REACTIVE OXIGEN", "REACTIVE OXYGEN", label, ignore.case = TRUE)
  ) %>%
  dplyr::arrange(desc(abs(NES))) %>%
  dplyr::mutate(label = factor(label, levels = rev(unique(label[order(NES)]))))

ggplot(df_mic, aes(x = label, y = NES, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = c("PBIC Enriched" = "#E66B6B", "MIC Enriched" = "#79C4D6")) +
  coord_flip() +
  labs(title = "Malignant PC PBIC vs MIC", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5))
