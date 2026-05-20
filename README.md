# SCdetMito

SCdetMito is an R package for sample-aware adaptive mitochondrial quality control in single-cell RNA-seq data.

Version: `1.3.6`  
Release date: `2026-05-20`

## Overview

SCdetMito provides mitochondrial QC decision support for Seurat workflows. It
constructs retained-cell profiles across candidate mitochondrial cutoffs,
quantifies interval-specific cell loss, estimates sample-level cutoffs using
retention-loss enrichment, and supports multi-sample QC strategies,
benchmarking, sensitivity analysis, cutoff confidence reporting, and provenance
tracking.

The package does not define a universal biological mitochondrial cutoff.
Detected cutoffs should be interpreted with retained-cell profiles,
interval-loss profiles, cutoff source, cutoff confidence, and study context.

## Installation

```r
install.packages("remotes")
remotes::install_github("ericHuu/SCdetMito@v1.3.6", build_vignettes = FALSE)
```

For local source installation:

```r
install.packages(
  "/path/to/SCdetMito_1.3.6.tar.gz",
  repos = NULL,
  type = "source"
)
```

## Example Data Options

SCdetMito includes a lightweight built-in demo object for quick checks. For
public real data, users can optionally load PBMC3K, the 10x Genomics mouse
heart 10k dataset, or the SeuratData IFNB-stimulated/control PBMC dataset.
PBMC3K is useful as a small smoke-test example, heart10k is the recommended
larger real-tissue example, and IFNB is useful for demonstrating
multi-condition QC strategies.

Online datasets are not included in the package. They are downloaded or loaded
only when explicitly called. Standard tests and package checks do not require
internet access.

## Built-In Lightweight Demo

`load_demo_seurat()` requires no internet access and is intended for smoke
tests, package examples, and quick output checks.

```r
library(SCdetMito)

seu <- load_demo_seurat()
det <- SCdetMito(
  seu,
  sample_col = "sample",
  mito_col = "mitoRatio",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

det$sample_cutoff_summary
```

## PBMC3K Public Smoke Test

`load_pbmc3k_online()` optionally loads the public 10x Genomics PBMC3K dataset.
PBMC3K is provided as a small public smoke-test example. Because PBMCs often
have relatively low mitochondrial fractions and dataset-specific QC structure,
PBMC3K cutoff results should be interpreted as a functionality check rather
than evidence of a universal mitochondrial threshold.

```r
pbmc <- load_pbmc3k_online(method = "10x")
det_pbmc <- SCdetMito(
  pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)
```

PBMC3K is downloaded only when `load_pbmc3k_online()` is called. Internet
access is required for `method = "10x"`. `SeuratData` is required for
`method = "SeuratData"` and the dataset must already be installed through
`SeuratData::InstallData("pbmc3k")`.

## Heart10k Real-Tissue Example

`load_heart10k_online()` optionally loads the public 10x Genomics mouse heart
10k dataset. It is the recommended larger real-tissue example for SoftwareX
demonstration of mitochondrial QC behavior in solid tissue. The dataset is
larger than the built-in demo and PBMC3K and may take more time and memory.

```r
heart <- load_heart10k_online()
det_heart <- SCdetMito(
  heart,
  sample_col = "sample",
  mito_col = "mitoRatio",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

det_heart$sample_cutoff_summary
```

## IFNB Multi-Condition Example

`load_ifnb_online()` optionally loads the public IFNB-stimulated/control PBMC
dataset from SeuratData. It is intended for demonstrating consensus, strictest,
and groupwise multi-sample/group QC strategies. IFNB requires SeuratData and
the dataset must already be installed through `SeuratData::InstallData("ifnb")`.

```r
ifnb <- load_ifnb_online()

qc_groupwise <- SCQCmulti(
  ifnb,
  sample_col = "sample",
  group_col = "condition",
  mito_col = "mitoRatio",
  cutoff_strategy = "groupwise",
  min_genes = 0,
  min_counts = 0,
  removeDouble = FALSE,
  write_plots = FALSE,
  write_tables = FALSE
)
```

## Main Functions

| Function | Purpose |
| --- | --- |
| `load_demo_seurat()` | Load a lightweight built-in demo object for smoke tests |
| `load_pbmc3k_demo()` | Load a small bundled PBMC3K-derived object for package examples |
| `load_pbmc3k_online()` | Optionally load PBMC3K for small public-data smoke testing |
| `load_heart10k_online()` | Optionally load the public 10x mouse heart 10k dataset for real-tissue examples |
| `load_ifnb_online()` | Optionally load IFNB-stimulated/control PBMC data for multi-condition QC examples |
| `add_mitoRatio()` | Add mitochondrial fraction to Seurat metadata |
| `SCdetMito()` | Estimate adaptive mitochondrial QC cutoffs |
| `SCQCone()` | Apply QC to one Seurat object |
| `SCQCmulti()` | Apply sample-aware multi-sample QC |
| `SCQCbenchmark()` | Compare fixed and adaptive QC strategies |
| `SCdetMito_sensitivity()` | Evaluate cutoff stability across parameter settings |
| `SCdetMito_methods()` | List available detection methods |

## Multi-Sample QC Strategies

- `consensus`: one sample-supported global cutoff for the integrated object.
- `strictest`: the minimum detected sample-level cutoff.
- `groupwise`: group-specific sample-supported cutoffs for studies where
  biological groups, batches, disease states, tissues, or species may differ in
  baseline mitochondrial profiles.

When `groupwise` is used and a group contains fewer than two samples,
`SCQCmulti()` warns that groupwise sample-supported cutoffs are less stable for
single-sample groups.

## Benchmarking and Sensitivity Analysis

`SCQCbenchmark()` compares fixed thresholds and adaptive strategies under
user-adjustable scoring profiles. The returned recommendation is a strategy
ranked under the selected scoring profile, not a universal QC rule.

```r
benchmark <- SCQCbenchmark(
  seurat_obj = seu,
  sample_col = "sample",
  group_col = "group",
  min_genes = 0,
  min_counts = 0,
  removeDouble = FALSE,
  run_downstream = FALSE,
  plot = FALSE,
  output_dir = tempdir()
)

benchmark$recommended_strategies
```

`SCdetMito_sensitivity()` evaluates cutoff stability across user-specified
parameter settings.

```r
sensitivity <- SCdetMito_sensitivity(
  seurat_obj = seu,
  sample_col = "sample",
  mito_col = "mitoRatio",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

sensitivity$sensitivity_summary
```

## Provenance Tracking

Major functions that modify a Seurat object store QC provenance in:

```r
seurat_obj@misc$SCdetMito_QC
```

The provenance record includes package metadata, timestamp, input and output
cell counts, retention rate, parameters, cutoff plan, cutoff source, cutoff
confidence, loss test, sample cutoff method, fallback events, optional doublet
filtering status, and session metadata.

DoubletFinder-based doublet removal is optional and not part of the core
retention-loss enrichment method.

## Interpretation Notes

- SCdetMito provides mitochondrial QC decision support.
- It does not define a universal biological mitochondrial cutoff.
- Detected cutoffs should be interpreted together with retained-cell profiles,
  interval-loss profiles, cutoff source, and cutoff confidence.
- Fallback-derived cutoffs are explicitly labeled and should be interpreted
  cautiously.
- Mitochondrial filtering should be combined with other QC metrics such as
  `nFeature`, `nCount`, doublet detection, ambient RNA assessment, cell type
  composition, and biological context.
- Benchmark scores are user-adjustable summaries and do not imply a
  universally applicable QC strategy.
- `atlas_benchmark_manifest()` is a reference-only planning manifest and does
  not validate dataset availability, download status, or annotation quality.

## Citation

```r
citation("SCdetMito")
```

Zenodo DOI: `[Zenodo DOI to be added]`

When reporting SCdetMito results, cite the package version, detection settings,
sample cutoff summary, and whether cutoffs were statistically detected or
fallback-derived.

## Contact

Silu Hu  
husilu0902@gmail.com  
Repository: https://github.com/ericHuu/SCdetMito
