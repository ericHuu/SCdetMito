# SCdetMito

`SCdetMito` is an R package for sample-aware adaptive mitochondrial quality
control in single-cell RNA-seq data, with a primary focus on Seurat objects.

Current package timestamp: `2026-05-18`; version `1.3.5`.

## What SCdetMito does

`SCdetMito` provides transparent QC decision support rather than a claim that
one universal mitochondrial cutoff exists. It estimates sample-level
mitochondrial ratio cutoffs from retained-cell profiles, aggregates those
cutoffs for multi-sample objects, and records QC provenance in Seurat object
metadata.

Main workflows:

- `add_mitoRatio()`: calculate and store mitochondrial ratio metadata.
- `SCdetMito()`: estimate sample-level cutoffs using retention-loss enrichment.
- `SCQCone()`: apply QC to a single sample or object.
- `SCQCmulti()`: apply integrated multi-sample QC with consensus, strictest, or
  groupwise sample-supported cutoffs.
- `SCQCbenchmark()`: compare fixed thresholds with adaptive strategies under
  user-adjustable scoring profiles.
- `SCdetMito_sensitivity()`: audit cutoff stability across parameter settings.

## Why fixed cutoffs can be problematic

A fixed mitochondrial cutoff can be adequate for a narrow single-sample
workflow. In integrated studies, baseline mitochondrial distributions may differ
by sample, batch, disease state, tissue, or species. A single fixed threshold
can therefore remove valid cells from one sample while retaining low-quality
cells in another. `SCdetMito` makes the cutoff decision auditable by returning
the retained-cell profile, interval-specific cell loss, candidate filters,
adjusted p-values, fallback status, and cutoff confidence.

## Core concept

For each sample `s` and candidate mitochondrial cutoff `c`, `SCdetMito`
calculates the retained-cell profile:

```text
R_s(c) = number of cells in sample s with mitoRatio < c
```

Adjacent cutoffs define interval-specific cell loss:

```text
L_s(c_j) = R_s(c_{j-1}) - R_s(c_j)
```

The workflow is:

```text
retained-cell profile
  -> interval-specific cell loss
  -> retention-loss enrichment test
  -> sample-level cutoff
  -> sample-supported global or groupwise cutoff
```

The default detector is `loss_test = "mad_zscore"` with
`p_adj_method = "BH"` and `sample_cutoff_method = "largest_drop"`.
`mad_zscore` is the primary robust FDR-compatible mode. `empirical_tail`
remains available as a non-parametric, highly discrete sensitivity mode.

## Installation

```r
# install.packages("remotes")
remotes::install_local("/path/to/SCdetMito")
```

## Quick Start

```r
library(SCdetMito)

seu <- load_demo_seurat()

qc_obj <- SCQCmulti(
  seurat_obj = seu,
  sample_col = "sample",
  group_col = "group",
  mode = "all",
  cutoff_strategy = "groupwise",
  max_mito = "SCdetMito",
  min_genes = 0,
  min_counts = 0,
  remove_doublets = FALSE,
  write_plots = FALSE,
  write_tables = FALSE
)

ncol(qc_obj)
qc_obj@misc$SCdetMito_QC
```

## Main Functions

```r
methods <- SCdetMito_methods()
methods$loss_tests
```

```r
det <- SCdetMito(
  seurat_obj = seu,
  sample_col = "sample",
  mito_col = "mitoRatio",
  loss_test = "mad_zscore",
  p_adj_method = "BH",
  sample_cutoff_method = "largest_drop",
  write_plots = FALSE,
  write_tables = FALSE,
  return_details = TRUE
)

det$sample_cutoff_summary
det$interval_loss_table
```

`change_points` is retained in returned objects only as a legacy alias of
`interval_loss_table`. The legacy CSV name `mito_change_point_results.csv` may
also be written for older workflows; new workflows should use
`SCdetMito_interval_loss_table.csv`.

## Multi-Sample QC Strategies

- `consensus`: sample-supported global cutoff for the integrated object. It
  selects the largest sample-level cutoff supported by at least the requested
  fraction of all samples, even when a group column is present.
- `strictest`: conservative global cutoff using the lowest detected
  sample-level cutoff.
- `groupwise`: group-specific sample-supported cutoff. This is useful when
  biological groups, batches, disease states, tissues, or species may differ in
  baseline mitochondrial profiles.

When `groupwise` is used and a group contains fewer than two samples,
`SCQCmulti()` warns that groupwise sample-supported cutoffs are less stable for
single-sample groups.

## Benchmarking

`SCQCbenchmark()` compares fixed thresholds and adaptive strategies. Composite
scores are user-adjustable decision-support summaries; they do not imply a
universal QC optimum.

Built-in scoring profiles:

- `balanced`
- `retention_focused`
- `stringent_mito`
- `balance_focused`

```r
benchmark <- SCQCbenchmark(
  seurat_obj = seu,
  sample_col = "sample",
  group_col = "group",
  min_genes = 0,
  min_counts = 0,
  remove_doublets = FALSE,
  run_downstream = FALSE,
  write_plots = FALSE,
  output_dir = tempdir()
)

benchmark$raw_benchmark_table
benchmark$normalized_score_table
benchmark$recommended_strategies
```

## Sensitivity Analysis

```r
sensitivity <- SCdetMito_sensitivity(
  seurat_obj = seu,
  sample_col = "sample",
  mito_col = "mitoRatio",
  bin_width = c(0.005, 0.01, 0.02),
  alpha = c(0.01, 0.05, 0.10),
  min_drop_fraction = c(0.005, 0.01, 0.02),
  loss_test = c("mad_zscore", "empirical_tail"),
  sample_cutoff_method = c("largest_drop", "max_significant"),
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

sensitivity$sensitivity_summary
```

The sensitivity summary reports sample-level cutoff median, IQR, min, max,
fallback fraction, not-detected fraction, and low-confidence fraction.

## Provenance Tracking

Major functions that modify a Seurat object store QC provenance in:

```r
seurat_obj@misc$SCdetMito_QC
```

The provenance record includes package metadata, timestamp, input/output cell
counts, retention rate, parameters, cutoff plan, cutoff source, cutoff
confidence, loss test, sample cutoff method, fallback events, optional
DoubletFinder status, and session metadata.

DoubletFinder-based doublet removal is optional and not part of the core
retention-loss enrichment method. If DoubletFinder is missing or fails, the
mitochondrial QC-filtered object is returned with a warning unless
`fail_action = "stop"`.

## Bundled Example Data

The package includes lightweight example objects for documentation and tests:

- `load_demo_seurat()`: synthetic multi-group Seurat object.
- `load_pbmc3k_demo()`: small PBMC3K-derived Seurat object for public-data
  smoke tests.

These data are intended for package examples, not for large-scale method
claims.

## Interpretation Notes

- `cutoff_confidence` is QC decision confidence, not biological certainty.
- Fallback-derived cutoffs are explicitly labeled and should be interpreted
  cautiously.
- `empirical_tail` p-values are rank-based and highly discrete; multiple
  testing adjustment may be conservative under dense cutoff scanning.
- `SCQCbenchmark()` ranks strategies under the selected scoring profile only.
- `atlas_benchmark_manifest()` is a reference-only planning manifest and does
  not validate dataset availability, download status, or annotation quality.

## Citation and Development Status

`SCdetMito` is under active validation for software-style submission and
reproducible review. Please cite the package version, selected settings, and
returned cutoff summaries when reporting results.
