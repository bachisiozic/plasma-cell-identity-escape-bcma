# `data/` — input tables

Place the input files here and update the paths at the top of each script (the originals
point to local Box / OneDrive locations). Files marked **controlled access** contain
patient-level data and should be obtained from the appropriate repository rather than
committed to GitHub.

## Genomic analysis

Used by `scripts/02_genomic_analysis.R` (analysis-ready tables produced by
`scripts/01_genomic_data_preparation.R`):

| File | Description | Sharing |
|---|---|---|
| `all_myeloma_and_cancer_driver_genes.txt` | Driver-gene coordinate reference (myeloma + cancer drivers, e.g. Rustad 2020) | shareable |
| `all_tce_cart_clin.txt` | Clinical/outcome table (PFS, OS, cohort, etc.) | controlled access |
| `jco_matrix_cart_tce.txt` | Maura et al. *JCO* 2024 genomic classification matrix | controlled access |
| `ALL_CNV_FINAL.txt` | Merged copy-number calls (all cohorts) | controlled access |
| `ALL_NS_SNV_FINAL.txt` | Non-synonymous mutations (dN/dS annotated) | controlled access |
| `ALL_SV_FINAL.txt` | Structural variants | controlled access |
| `final_matrix_genomic_drivers.txt` | Per-sample driver matrix | controlled access |
| `biallelic_heatmap_2025.txt` | Biallelic BCMA-region loss annotation (*TRAF3/TRAF2/CYLD*) | controlled access |
| `CAR-T_RNAseq_Data.xlsx` or `rna_seq_bcma.txt` | Bulk RNA-seq (z-scored *TNFRSF17/TP53/CD38/XBP1/SP140*) | shareable (de-identified) |

`scripts/01_genomic_data_preparation.R` additionally reads several raw per-cohort files
(`all_cnv_*`, `all_sv_*`, `MSK_TCE_*`, `output_dnds_all.RDS`). These are upstream and
**controlled access**; most users should start from the `ALL_*_FINAL.txt` tables above.

## Single-cell analysis

Used by `scripts/03_scRNAseq_analysis.Rmd` (large objects — deposit on Zenodo/GEO):

| File | Description |
|---|---|
| `Pt_Malignant_plasma_pt1_excluded_Pre_Tx.RDS` | Moffitt malignant plasma-cell Seurat object |
| `baseline_iTME.rds` | Moffitt immune-microenvironment Seurat object |
| `MM_CART_BCD_integrated50pc_..._20231119.rds` | Dhodapkar (public) integrated Seurat object |
| `number_of_cells_per_each_sample.txt` | Per-sample cell-type counts (immune composition) |

Gene sets for the GSEA step live in [`../gene_sets/`](../gene_sets/README.md).
