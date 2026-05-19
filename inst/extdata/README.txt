Bundled example data
====================

`demo_multigroup_seurat.rds` is a small toy Seurat object bundled only for
examples, vignette demonstration, and fast unit tests.

`pbmc3k_lightweight_seurat.rds` is a small subset derived from the public 10x
Genomics PBMC3K dataset. It is included to let users run SCdetMito examples on
a recognizable public dataset without shipping the full PBMC3K matrix in the
installed package.

It is not intended to replace validation on real public data. Real-data
validation scripts are stored under `data-raw/` and optional local validation
artifacts are stored under `data-public/`.
