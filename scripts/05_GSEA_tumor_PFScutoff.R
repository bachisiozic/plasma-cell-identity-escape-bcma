library(tidyr)
library(dplyr)
library(pheatmap)
library(Seurat)
library(edgeR)
library(DESeq2)

tumor2=readRDS(".../moffitt_collapsed_counts_only_plasmacells.RDS")
all.2=readRDS(".../dhodapkar_collapsed_counts_only_plasmacells.RDS")

#### . . . - Differential expression analysis
moffitt.exp=tumor2$pseudobulk_sum
all.exp=all.2$pseudobulk_sum
### . . - Selecting the genes that are shared by both cohorts
GENES=intersect(rownames(all.exp),rownames(moffitt.exp))

#### . . . - Creating the dataframe for moffitt with selected genes
moffitt.exp=tumor2$pseudobulk_sum

moffitt.exp=moffitt.exp[rownames(moffitt.exp) %in% GENES,]
moffitt.exp.t=t(moffitt.exp)

### . . - Excluding pt 11 because died in 30 days after infusion
moffitt.deg=merge(moffitt.df.2,moffitt.exp.t,by.x="sample",by.y="row.names",all.x=T)

#### . . . - Creating the dataframe for Dhodapkar with selected genes
all.exp=all.2$pseudobulk_sum
all.exp=all.exp[rownames(all.exp) %in% GENES,]
all.exp.t=t(all.exp)
all.deg=merge(all.df.1,all.exp.t,by.x="sample",by.y="row.names",all.x=T)

#### . . . - Creating one big dataframe integrating both cohorts
df.deg=plyr::rbind.fill(moffitt.deg,all.deg)

#### . . . - Selecting only the baseline samples
df.deg.bas=df.deg[df.deg$time_point %in% c("PreTx","Pre"),]

#### . . . Creatig input data
expr=df.deg.bas[,c(1,6:ncol(df.deg.bas))]
rownames(expr)=expr$sample
expr=expr[,-1]
expr=t(expr)

pdata=df.deg.bas[,1:5]
rownames(pdata)=pdata$sample
pdata=pdata[,-1]

info=read.delim(".../PFS_different_timepoints_Moffitt_Dhodapkar.txt",
                stringsAsFactors = F)
info=info[,-2]

pdata.1=merge(pdata,info,by.x="row.names",by.y="sample",all.x=T)
rownames(pdata.1)=pdata.1$Row.names

#### . . . - DEG analysis
NEW=NULL
for(j in c("D60_response","D90_response","D120_response","D150_response","D180_response","D270_response","D360_response")){
  dds <- DESeqDataSetFromMatrix(countData = expr,
                                colData = pdata.1,
                                design = as.formula(paste0("~ ", j)))
  
  keep <- rowSums(counts(dds) >= 10) >= ceiling(0.2*ncol(dds))
  dds  <- dds[keep,]
  dds  <- DESeq(dds)
  res  <- results(dds, contrast = c(j,"NDR","DR"))
  
  ranks <- res$stat 
  ok <- !is.na(ranks)
  ranks <- ranks[ok]
  ranks[ranks == Inf] <- 999
  names(ranks) <- rownames(res)[ok]
  
  
  CIN70 <- c("TPX2", "PRC1", "FOXM1", "CDC2", "C20orf24", "TGIF2", "MCM2", "H2AFZ",
             "TOP2A", "PCNA", "UBE2C", "MELK", "TRIP13", "CNAP1", "MCM7", "RNASEH2A",
             "RAD51AP1", "KIF20A", "CDC45L", "MAD2L1", "ESPL1", "CCNB2", "FEN1", "TTK",
             "CCT5", "RFC4", "ATAD2", "ch-TOG", "NUP205", "CDC20", "CKS2", "RRM2", "ELAVL1",
             "CCNB1", "RRM1", "AURKB", "MSH6", "EZH2", "CTPS", "DKC1", "OIP5", "CDCA8",
             "PTTG1", "CEP55", "H2AFX", "CMAS", "BRRN1", "MCM10", "LSM4", "MTB", "ASF1B",
             "ZWINT", "TOPK", "FLJ10036", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2", "RAD21",
             "ACTL6A", "GPI", "MGC13096", "SFRS2", "HDGF", "NXT1", "NEK2", "DHCR7", "STK6",
             "Aurora-A", "NDUFAB1", "KIAA0286", "KIF4A")
  
  library(fgsea)
  library(msigdbr)
  setlist<-read.table(".../h.all.v6.2.symbols.gmt_MOD_CART_paper.txt",sep="\t",fill=T,header=F,as.is=T)
  gs<-vector("list",length=nrow(setlist))
  names(gs)<-setlist[,1]
  setlist<-setlist[,-c(1:2)]
  for(i in 1:nrow(setlist)){
    temp<-as.vector(as.matrix(setlist[i,]))
    temp<-temp[temp!=""]
    gs[[i]]<-temp
  }
  
  bm.pc=read.delim(".../HAY_BONE_MARROW_PLASMA_CELL.v2025.1.Hs.gmt",
                   sep="\t",fill=T,header=F,as.is=T)
  gs1<-vector("list",length=nrow(bm.pc))
  names(gs1)<-bm.pc[,1]
  bm.pc<-bm.pc[,-c(1:2)]
  for(i in 1:nrow(bm.pc)){
    temp<-as.vector(as.matrix(bm.pc[i,]))
    temp<-temp[temp!=""]
    gs1[[i]]<-temp
  }
  
  sets  <- c(gs, list(CIN70=CIN70),gs1)
  fg <- fgsea(pathways = sets, stats = ranks,nperm = 10000)
  NDR=fg[order(fg$padj), c("pathway","NES","padj")]
  
  write.table(as.data.frame(NDR),
              paste(".../",
                    j,
                    "_Malignant_sc_plasma_cells_GSEA.txt",
                    sep=""),
              col.names = T,
              row.names = F,
              quote = F,
              sep="\t")
  
  NEW=rbind(NEW,
            data.frame(pathway=NDR$pathway,
                       NES=NDR$NES,
                       FDR=NDR$padj,
                       PFS=j))
}

NEW.1=NEW

# --- 1. Enforce the correct column (time-point) order --------------------
pfs_levels <- c("D60_response", "D90_response",
                "D120_response", "D150_response", "D180_response",
                "D270_response", "D360_response")
NEW.1$PFS <- factor(NEW.1$PFS, levels = pfs_levels)

# --- 2. Long -> wide: NES matrix (pathways = rows, PFS = cols) ------------
nes_mat <- NEW.1 %>%
  select(pathway, PFS, NES) %>%
  pivot_wider(names_from = PFS, values_from = NES) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# Reorder columns to the intended time order (pivot_wider follows factor
# levels, but this makes it explicit and safe)
nes_mat <- nes_mat[, pfs_levels]

# --- 3. Same reshape for FDR, then build the asterisk matrix -------------
fdr_mat <- NEW.1 %>%
  select(pathway, PFS, FDR) %>%
  pivot_wider(names_from = PFS, values_from = FDR) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()
fdr_mat <- fdr_mat[rownames(nes_mat), pfs_levels]   # align to nes_mat exactly
star_mat <- ifelse(!is.na(fdr_mat) & fdr_mat < 0.1, "*", "")

# --- 4. Symmetric diverging color scale (blue < 0 < red, white at 0) -----
max_abs <- max(abs(nes_mat), na.rm = TRUE)
breaks  <- seq(-max_abs, max_abs, length.out = 101)   # 100 colors -> 101 breaks
palette <- colorRampPalette(c("#79C4D6", "white", "#E66B6B"))(100)

# --- 5. Draw the heatmap -------------------------------------------------
pheatmap(nes_mat,
         color            = palette,
         breaks           = breaks,
         cluster_rows     = TRUE,    
         cluster_cols     = FALSE,   
         display_numbers  = star_mat,
         number_color     = "black",
         fontsize_number  = 12,
         na_col           = "grey90",
         border_color     = NA,
         main             = "NES across time points (* FDR < 0.1)")


NEW.1$n.NDR.patients=ifelse(NEW.1$PFS %in% c("D60_response","D90_response"),5,
                   ifelse(NEW.1$PFS %in% c("D120_response","D150_response","D180_response","D270_response"),6,
                          ifelse(NEW.1$PFS %in% c("D360_response"),7,0)))

write.table(NEW.1,
            ".../GSEA_different_PFS_cutoff.txt",
            col.names = T,
            row.names = F,
            quote = F,
            sep="\t")

df.bas.1=merge(df.bas,info,by="sample",all.x=T)

BCMA=NULL
for(i in c(49,50,51,52,53,54,4,55)){
  test=wilcox.test(df.bas.1$TNFRSF17 ~ df.bas.1[,i],alternative = "greater")
  BCMA=rbind(BCMA,
              data_frame(cut.off=colnames(df.bas.1[i]),
                         p.value=test$p.value))
}

# A tibble: 8 × 2
#   cut.off       p.value
#   <chr>           <dbl>
# 1 D60_response  0.00962
# 2 D90_response  0.00962
# 3 D120_response 0.0363
# 4 D150_response 0.0363
# 5 D180_response 0.0363
# 6 D270_response 0.0363
# 7 D360_response 0.0571



