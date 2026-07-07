# Plasma cell identity escape drives resistance to anti-BCMA T-cell–redirecting therapy in multiple myeloma

- Bachisio Ziccheddu — `ziccheb@mskcc.org`
- Francesco Maura — `mauraf@mskcc.org`

This repository contains the analysis code accompanying the manuscript:

> **Plasma cell identity escape drives resistance to anti-BCMA T-cell–redirecting therapy in multiple myeloma.**
> *Maura F., et al.

The study integrates three complementary data modalities to understand why some multiple
myeloma patients lose durable response to BCMA-targeted T-cell–redirecting therapy
(CAR-T and bispecific antibodies, "TCE"):

1. **Whole-genome sequencing (WGS) genomics** — copy-number, structural variants,
   non-synonymous mutations and driver genes; survival modelling; biallelic *TNFRSF17*
   (BCMA) loss; and a combined "complex genomics + plasma-cell-identity loss" signature.
2. **Single-cell RNA-seq (scRNA-seq)** — the malignant plasma-cell compartment and the
   surrounding immune microenvironment, comparing durable (DR) versus non-durable (NDR)
   responders.
3. **Mouse model validation** — two independent mouse models (Spain cohort: BIC/PBIC/MIC
   samples; vkMYC model: responders vs non-responders) confirming plasma-cell identity
   loss and GSEA signatures observed in human data.

Throughout, **DR = durable response** (progression-free survival > 180 days) and
**NDR = non-durable response**.

---

## Repository structure

```
.
├── README.md
├── LICENSE
├── .gitignore
├── scripts/
│   ├── 01_genomic_analysis.R           # WGS analyses + main/supplementary figures
│   ├── 02_scRNAseq_analysis.Rmd        # scRNA-seq report (tumor + immune), knit to HTML
│   ├── 03_spain_mouse_analysis.R       # Spain cohort mouse scRNA-seq + GSEA (BIC/PBIC/MIC)
│   └── 04_vkmyc_mouse_analysis.R       # vkMYC mouse model DEG analysis + GSEA
│   └── 05_GSEA_tumor_PFScutoff.R       # GSEA using different PFS thresholds
│   └── 06_immune_cells_PFScutoff.R     # Looking the T cells compartment using different PFS thresholds
└── gene_sets/                          # gene sets used for GSEA / enrichment
   ├── README.md
   └── CIN70.gmt
   └── h.all.v6.2.symbols.gmt_MOD_CART_paper.txt
   └── HAY_BONE_MARROW_PLASMA_CELL.v2025.1.Hs.gmt
   └── Mus_musculus_hallmarks_cin70_plasmacell.RData
```

---

## Software and dependencies

Analyses were run in **R (≥ 4.2)**. Install the packages below before running.

**Genomic analysis (`01`–`02`)**

```r
install.packages(c("stringr","dplyr","plyr","tidyr","tidyverse","data.table",
                   "ggplot2","survival","survminer","RColorBrewer","colorspace",
                   "reshape","readxl","pheatmap","forestmodel","deconstructSigs",
                   "quantsmooth"))
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges","IRanges","GenomicFeatures","maftools",
                       "ComplexHeatmap","biomaRt"))
# specialised (myeloma) packages
# mmsig:  https://github.com/bachisiozic/mmsig
# dndscv:  https://github.com/im3sanger/dndscv
```

**Single-cell analysis (`03`)**

```r
install.packages(c("Seurat","SeuratObject","harmony","dplyr","stringr","tidyr","ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","edgeR","fgsea","msigdbr"))
```

**Mouse model analyses (`04`–`05`)**

```r
install.packages(c("Seurat","dplyr","ggplot2","patchwork","ggpubr","rstatix","reshape2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","edgeR","fgsea","msigdbr","ComplexHeatmap",
                       "clusterProfiler","EnhancedVolcano"))
```

---



## Data availability

- WGS: `EGAS50000001817` and `EGAS50000000546`

---
