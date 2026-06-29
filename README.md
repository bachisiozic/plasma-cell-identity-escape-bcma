# Plasma cell identity escape drives resistance to anti-BCMA T-cell–redirecting therapy in multiple myeloma

This repository contains the analysis code accompanying the manuscript:

> **Plasma cell identity escape drives resistance to anti-BCMA T-cell–redirecting therapy in multiple myeloma.**
> *Maura F., et al.

The study integrates two complementary data modalities to understand why some multiple
myeloma patients lose durable response to BCMA-targeted T-cell–redirecting therapy
(CAR-T and bispecific antibodies, "TCE"):

1. **Whole-genome sequencing (WGS) genomics** — copy-number, structural variants,
   non-synonymous mutations and driver genes; survival modelling; biallelic *TNFRSF17*
   (BCMA) loss; and a combined "complex genomics + plasma-cell-identity loss" signature.
2. **Single-cell RNA-seq (scRNA-seq)** — the malignant plasma-cell compartment and the
   surrounding immune microenvironment, comparing durable (DR) versus non-durable (NDR)
   responders.

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
│   ├── 01_genomic_data_preparation.R   # build the analysis-ready genomic tables
│   ├── 02_genomic_analysis.R           # WGS analyses + main/supplementary figures
│   └── 03_scRNAseq_analysis.Rmd        # scRNA-seq report (tumor + immune), knit to HTML
├── gene_sets/                          # gene sets used for GSEA / enrichment
│   ├── README.md
│   └── CIN70.gmt
│   └── h.all.v6.2.symbols.gmt_MOD_CART_paper.txt
│   └── HAY_BONE_MARROW_PLASMA_CELL.v2025.1.Hs.gmt
│   └── Mus_musculus_hallmarks_cin70_plasmacell.RData
├── data/                               # input tables (see data/README.md)
│   └── README.md
└── results/                            # derived "source data" tables that reproduce figures
    └── README.md
```

The two modalities are **independent**: you can run the genomic analysis without the
single-cell analysis and vice versa.

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
# mmsig:   https://github.com/UM-Myeloma-Genomics/mmSig
# dndscv:  https://github.com/im3sanger/dndscv
```

**Single-cell analysis (`03`)**

```r
install.packages(c("Seurat","SeuratObject","harmony","dplyr","stringr","tidyr","ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","edgeR","fgsea","msigdbr"))
```

---


## Data availability

- WGS / variant calls: `<EGA/dbGaP accession>`

---


## Contact

- Bachisio Ziccheddu — `ziccheb@mskcc.org`
- Francesco Maura — `mauraf@mskcc.org`

