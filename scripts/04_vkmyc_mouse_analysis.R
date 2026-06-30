# =============================================================================
# vkMYC mouse model — DEG analysis and GSEA
# =============================================================================
# Analysis outline:
#   Part 1 — Data loading and sample-level summarization
#   Part 2 — Pseudobulk differential expression (DESeq2)
#   Part 3 — GSEA using MSigDB Hallmarks + CIN70 + plasma-cell identity gene set
#   Part 4 — Visualization
# =============================================================================

library(edgeR)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(fgsea)
library(msigdbr)


# =============================================================================
# PART 1 — Load data and collapse technical replicates per sample
# =============================================================================

rna <- read.delim(
  ".../RNA_fra_annot.txt",
  stringsAsFactors = FALSE
)

# Remove columns labelled "V_V" (variant-to-variant technical duplicates)
rna1 <- rna[, grep("V_V", colnames(rna), invert = TRUE)]

# Sample metadata: 10 vkMYC mice with known treatment outcome
sample_info <- data.frame(
  sample  = c("Vk40149","Vk39225","Vk32908","Vk12653","Vk12598",
               "Vk29790","Vk39814","Vk27181","Vk35633","Vk35984"),
  outcome = c(rep("resp", 7), rep("nonresp", 3))
)
dim(sample_info)  # 10 x 2

# Collapse multiple columns per sample into a single median column
# (some samples have >1 technical replicate in the expression matrix)
for (s in sample_info$sample) {
  matches <- colnames(rna1)[grep(s, colnames(rna1))]
  if (length(matches) == 0) next

  rna1[, paste0(s, ".median")] <- if (length(matches) > 1) {
    apply(rna1[, matches, drop = FALSE], 1, median, na.rm = TRUE)
  } else {
    rna1[, matches]
  }
}

# Keep gene annotation columns (1:4) and the newly created median columns
rna2 <- rna1[, c(1:4, grep("median", colnames(rna1)))]
dim(rna2)

# Retain SYMBOL column and median expression; drop unannotated rows
expr <- rna2[, c(2, 5:10)]
expr <- expr[complete.cases(expr$SYMBOL), ]


# =============================================================================
# PART 2 — Pseudobulk differential expression: non-responders vs responders
# =============================================================================

# Convert to matrix and clean column names
rna.mat <- as.matrix(expr[, -1])
rownames(rna.mat) <- expr$SYMBOL
colnames(rna.mat) <- gsub(".median", "", colnames(rna.mat))

# Round to integers for DESeq2 (values are counts or count-like)
expr.counts <- round(rna.mat)

# Align sample metadata to matrix columns
pdata <- sample_info[sample_info$sample %in% colnames(expr.counts), ]
rownames(pdata) <- pdata$sample
pdata <- pdata[, -1, drop = FALSE]

# DESeq2 pipeline
dds  <- DESeqDataSetFromMatrix(countData = expr.counts, colData = pdata, design = ~ outcome)
keep <- rowSums(counts(dds) >= 10) >= ceiling(0.2 * ncol(dds))
dds  <- dds[keep, ]
dds  <- DESeq(dds)

# Contrast: non-responders vs responders
res <- results(dds, contrast = c("outcome","nonresp","resp"))

# Spot-check key myeloma genes
res[rownames(res) == "Tnfrsf17", ]
res[rownames(res) == "Xbp1", ]
res[rownames(res) == "Myc", ]

# Save full DEG table
write.table(as.data.frame(res),
  ".../DEGs_nonresp_VS_resp_vkMYC.txt",
  col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")


# =============================================================================
# PART 3 — GSEA
# =============================================================================

# Load pre-computed DEG results (skip if running sequentially from Part 2)
res <- read.delim(
  ".../DEGs_nonresp_VS_resp_vkMYC.txt",
  stringsAsFactors = FALSE
)

# Rank genes by log2 fold change (non-responders vs responders)
ranks <- res$stat
ok    <- !is.na(ranks)
ranks <- ranks[ok]
ranks[ranks == Inf] <- 999
names(ranks) <- rownames(res)[ok]

# --- Gene sets ---------------------------------------------------------------
# CIN70: chromosomal instability signature (mouse gene symbols)
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

# Plasma-cell identity gene set (mouse orthologs of known PC markers)
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

# MSigDB Hallmark gene sets (mouse)
m_df     <- msigdbr(species = "Mus musculus", category = "H")
pathways <- split(m_df$gene_symbol, m_df$gs_name)

sets <- c(pathways,
          list(CIN70       = CIN70_mouse,
               PLASMA_CELL = plasma_cell))

# Run GSEA
fg       <- fgsea(pathways = sets, stats = ranks, nperm = 10000)
res_gsea <- fg[order(fg$padj), c("pathway","NES","padj")]


# =============================================================================
# PART 4 — Visualization
# =============================================================================

# Enrichment plot for the plasma-cell identity gene set
plotEnrichment(pathway = sets[["PLASMA_CELL"]], stats = ranks) +
  labs(title = "PLASMA CELL Identity", subtitle = "Non-responders vs Responders") +
  theme_minimal()

# Bar plot of significant pathways (padj < 0.1)
sig_gsea <- as.data.frame(res_gsea[res_gsea$padj < 0.1, ])

df_plot <- sig_gsea %>%
  dplyr::mutate(
    direction = ifelse(NES > 0, "Non-responder Enriched", "Responder Enriched"),
    label     = gsub("^HALLMARK_", "", pathway),
    label     = gsub("_", " ", label),
    label     = gsub("REACTIVE OXIGEN", "REACTIVE OXYGEN", label, ignore.case = TRUE)
  ) %>%
  dplyr::arrange(desc(abs(NES))) %>%
  dplyr::mutate(label = factor(label, levels = rev(unique(label[order(NES)]))))

ggplot(df_plot, aes(x = label, y = NES, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = c("Non-responder Enriched" = "#E66B6B",
                               "Responder Enriched"     = "#79C4D6")) +
  coord_flip() +
  labs(title = "vkMYC — Non-responders vs Responders", x = NULL, y = "NES") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.title       = element_blank(),
    plot.title         = element_text(face = "bold", hjust = 0.5)
  )
