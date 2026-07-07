library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(Seurat)

#### . . . - Importing the data
baseline1=readRDS(".../baseline.1.RDS")
info=read.delim(".../PFS_different_timepoints_Moffitt_Dhodapkar.txt",
                stringsAsFactors = F)
info=info[,-1]
colnames(info)[1]="sample"
info=info[,-2]


rownames(info) <- info$sample
cell_samples <- baseline1$sample
new_meta <- info[cell_samples, setdiff(colnames(info), "sample")]
rownames(new_meta) <- colnames(baseline1)   # cell barcodes, in order

# Add it to the object
baseline1 <- AddMetaData(baseline1, metadata = new_meta)

#### . . . - CD8 T cells
for(j in c("D60_response","D90_response","D120_response","D150_response","D180_response","D270_response","D360_response")){
  Idents(baseline1) <- "predicted.celltype.l1"   
  de_plasma <- FindMarkers(
    baseline1,
    ident.1 = "NDR", ident.2 = "DR",
    group.by = j,   # compare DR vs NDR
    subset.ident = c("CD8 T"),
    test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1
  )
  
  new=de_plasma[rownames(de_plasma) %in% c("CD160", 
                                           "CST7",
                                           "GZMB", "HAVCR2",
                                           "LAG3",
                                           "CTLA4","TCF7","EOMES", 
                                           "PRDM1","BATF","IFNG","TNF","HLA-DR","GITR","CD57","IRF4", 
                                           "BTLA","ICOS", 
                                           "4-1BB","PRF1","CD38", 
                                           "KLRG1", 
                                           "PDCD1","TIGIT",
                                           "KLRD1",
                                           "PTPN12","HLA-DRB1",
                                           "TOX","IL7R","CD69","FOS"),]
  
  new$gene=rownames(new)
  
  de_plasma$gene=rownames(de_plasma)
  
  pstars <- function(p) {
    ifelse(p < 1e-4, "****",
           ifelse(p < 1e-3, "***",
                  ifelse(p < 1e-2, "**",
                         ifelse(p < 0.05, "*", ""))))
  }
  
  df_long <- new %>%
    transmute(
      gene,
      FDR = p_val_adj,
      lfc = avg_log2FC,
      NDR_pct = pct.1 * 100,
      DR_pct  = pct.2 * 100
    ) %>%
    pivot_longer(c(NDR_pct, DR_pct), names_to = "group", values_to = "pct") %>%
    mutate(
      cutoff = ifelse(group == "NDR_pct", "NDR", "DR"),
      # color encodes enrichment direction (mirror sign for DR so red/blue match the example)
      color_val = ifelse(cutoff == "NDR", lfc, -lfc),
      stars = pstars(FDR),
      p_label = sprintf("FDR=%.2g", FDR)
    )
  
  
  gene_order <- rev(df_long %>% distinct(gene) %>% pull(gene))
  df_long$gene <- factor(df_long$gene, levels = gene_order)
  
  ggplot(df_long, aes(x = cutoff, y = gene)) +
    geom_point(aes(size = pct, fill = color_val),
               shape = 21, colour = "grey25", stroke = 0.25) +
    # add significance stars above the dots
    geom_text(aes(label = stars),
              vjust = -0.9, size = 3.5, lineheight = 0.9, na.rm = TRUE) +
    scale_size_continuous("Percent Expressed", range = c(1.8, 7),
                          breaks = c(25, 50, 75), limits = c(0, 100)) +
    scale_fill_gradient2("log2(FC)",
                         low = "#2C7FB8", mid = "white", high = "#D7301F",
                         midpoint = 0, limits = c(-4, 4)) +
    labs(title = paste("Genes related to exhaustion - CD8 T - ", j,sep=""), x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}


#### . . . - CD4 T cells
for(j in c("D60_response","D90_response","D120_response","D150_response","D180_response","D270_response","D360_response")){
  Idents(baseline1) <- "predicted.celltype.l1"
  de_plasma <- FindMarkers(
    baseline1,
    ident.1 = "NDR", ident.2 = "DR",
    group.by = j,   # compare DR vs NDR
    subset.ident = c("CD4 T"),
    test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1
  )
  
  new=de_plasma[rownames(de_plasma) %in% c("CD160", 
                                           "CST7",
                                           "GZMB", "HAVCR2",
                                           "LAG3",
                                           "CTLA4","TCF7","EOMES", 
                                           "PRDM1","BATF","IFNG","TNF","HLA-DR","GITR","CD57","IRF4", 
                                           "BTLA","ICOS", 
                                           "4-1BB","PRF1","CD38", 
                                           "KLRG1", 
                                           "PDCD1","TIGIT",
                                           "KLRD1",
                                           "PTPN12","HLA-DRB1",
                                           "TOX","IL7R","CD69","FOS"),]
  
  
  
  new$gene=rownames(new)
  
  de_plasma$gene=rownames(de_plasma)
  
  pstars <- function(p) {
    ifelse(p < 1e-4, "****",
           ifelse(p < 1e-3, "***",
                  ifelse(p < 1e-2, "**",
                         ifelse(p < 0.05, "*", ""))))
  }
  
  df_long <- new %>%
    transmute(
      gene,
      FDR = p_val_adj,
      lfc = avg_log2FC,
      NDR_pct = pct.1 * 100,
      DR_pct  = pct.2 * 100
    ) %>%
    pivot_longer(c(NDR_pct, DR_pct), names_to = "group", values_to = "pct") %>%
    mutate(
      cutoff = ifelse(group == "NDR_pct", "NDR", "DR"),
      # color encodes enrichment direction (mirror sign for DR so red/blue match the example)
      color_val = ifelse(cutoff == "NDR", lfc, -lfc),
      stars = pstars(FDR),
      p_label = sprintf("FDR=%.2g", FDR)
    )
  
  
  gene_order <- rev(df_long %>% distinct(gene) %>% pull(gene))
  df_long$gene <- factor(df_long$gene, levels = gene_order)
  
  ggplot(df_long, aes(x = cutoff, y = gene)) +
    geom_point(aes(size = pct, fill = color_val),
               shape = 21, colour = "grey25", stroke = 0.25) +
    # add significance stars above the dots
    geom_text(aes(label = stars),
              vjust = -0.9, size = 3.5, lineheight = 0.9, na.rm = TRUE) +
    scale_size_continuous("Percent Expressed", range = c(1.8, 7),
                          breaks = c(25, 50, 75), limits = c(0, 100)) +
    scale_fill_gradient2("log2(FC)",
                         low = "#2C7FB8", mid = "white", high = "#D7301F",
                         midpoint = 0, limits = c(-4, 4)) +
    labs(title = paste("Genes related to exhaustion - CD4 T - ", j,sep=""), x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}







