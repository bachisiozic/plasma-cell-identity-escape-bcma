# `gene_sets/` — gene sets for GSEA / enrichment

The single-cell GSEA step (`scripts/03_scRNAseq_analysis.Rmd`, Part 1) uses the gene sets
below. Place the `.gmt` / `.txt` files here and point the `read.table()` / `read.delim()`
calls in the report to this folder.

| File | Description | Source |
|---|---|---|
| `CIN70.gmt` | 70-gene chromosomal-instability / proliferation signature | provided here (defined in the analysis code) |
| `h.all.symbols.gmt` | MSigDB **Hallmark** collection (paper-specific subset used in the manuscript) | MSigDB — https://www.gsea-msigdb.org/gsea/msigdb |
| `HAY_BONE_MARROW_PLASMA_CELL.gmt` | Bone-marrow plasma-cell identity signature | MSigDB (Hay et al.) — https://www.gsea-msigdb.org/gsea/msigdb |
| `IFNG_GS`, `ISG_RS` | Interferon-γ / interferon-stimulated-gene signatures | as cited in the manuscript |

## A note on redistribution

`CIN70.gmt` is included because the gene list is defined directly in the analysis code.
The **MSigDB**-derived sets (Hallmark, `HAY_BONE_MARROW_PLASMA_CELL`) are **not
redistributed here** — MSigDB has its own license/registration. Download them from MSigDB
and place them in this folder. The report expects plain `.gmt`-style tab-delimited files
(set name in column 1, an optional description in column 2, then one gene symbol per
column).

> **Tip for reproducibility:** record the exact MSigDB version you used (e.g. `v6.2`,
> `v2025.1.Hs`) in the manuscript methods, since gene-set contents change between releases.
