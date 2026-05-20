Bundled example data
====================

`demo_multigroup_seurat.rds` is a lightweight demonstration Seurat object
generated for package examples, vignette demonstration, and fast unit tests. It
should not be used for biological inference.

`pbmc3k_lightweight_seurat.rds` is a small subset derived from the public 10x
Genomics PBMC3K dataset. It is included for public-data smoke tests without
shipping the full PBMC3K matrix in the installed package. The object was
downsampled/processed for fast testing. Users requiring the full public PBMC3K
data should use `load_pbmc3k_online()`.

It is not intended to replace validation on larger public datasets. Optional
online loaders for PBMC3K, 10x mouse brain 5k, 10x mouse heart 10k, and
SeuratData IFNB are available as user-run examples.

No private, unpublished, or patient-identifiable data are included in these
bundled example objects.
