# SCdetMito

SCdetMito is an R package for reference-aware and sample-aware adaptive mitochondrial quality control in single-cell RNA-seq data.

Version: `1.4.0`  
Release date: `2026-05-20`

## Overview

SCdetMito provides mitochondrial QC decision support for Seurat workflows. It
constructs retained-cell profiles across candidate mitochondrial cutoffs,
quantifies interval-specific cell loss, estimates multiple interpretable
sample-level cutoff candidates, integrates literature-informed reference
cutoffs, and supports multi-sample QC strategies, benchmarking, sensitivity
analysis, cutoff confidence reporting, warning fields, and provenance tracking.

The package does not define a universal biological mitochondrial cutoff.
Detected cutoffs should be interpreted with retained-cell profiles,
interval-loss profiles, reference values, cutoff source, cutoff confidence, and
study context.

## Installation

```r
install.packages("remotes")
remotes::install_github("ericHuu/SCdetMito@v1.4.0", build_vignettes = FALSE)
```

For local source installation:

```r
install.packages(
  "/path/to/SCdetMito_1.4.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

## Reference-Aware Cutoff Interpretation

SCdetMito v1.4.0 separates several cutoff concepts that can differ in real
tissue datasets:

| Cutoff type | Meaning |
| --- | --- |
| `reference_cutoff` | Literature-informed species/tissue reference used as a decision prior |
| `first_significant_cutoff_high` | Most permissive significant retention-loss boundary when scanning from high to low cutoff |
| `largest_drop_cutoff` | Significant cutoff with the largest interval-specific cell loss |
| `selected_cutoff` | Cutoff selected under `sample_cutoff_method` |
| `detected_cutoff` | Backward-compatible alias of `selected_cutoff` |

`reference_cutoff` values are not hard rules. They help flag cases where a
data-driven cutoff is far above common literature expectations or where a
largest-drop cutoff is very stringent and retains few cells. When no exact
tissue-specific value is available, SCdetMito falls back to a species-level
reference and reports this explicitly; tissue-note rows should not be treated
as exact tissue-specific cutoffs unless `is_tissue_specific = TRUE`.

```r
SCdetMito_reference_cutoffs()

det <- SCdetMito(
  seu,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "human",
  tissue = "global",
  sample_cutoff_method = "reference_guided",
  return_details = TRUE
)

det$sample_cutoff_summary[, c(
  "sample",
  "reference_cutoff",
  "first_significant_cutoff_high",
  "largest_drop_cutoff",
  "selected_cutoff",
  "detected_cutoff",
  "warning_level",
  "has_reference_deviation"
)]
```

## Example Data Options

SCdetMito includes a lightweight built-in demo object for quick checks.
Optional public datasets can be loaded for real-data demonstrations. PBMC3K is
a small public smoke-test example. Mouse brain/nuclei and heart10k are
real-tissue candidate datasets for inspecting mitochondrial cutoff behavior.
Heart10k may have high mitochondrial burden and should be interpreted as a
stress-test. IFNB is a multi-condition PBMC dataset useful for demonstrating
multi-sample/group QC strategies.

Online datasets are not included in the package. They are downloaded or loaded
only when explicitly called. Standard tests and package checks do not require
internet access. Public dataset suitability should be checked before using
figures in a manuscript.

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
  species = "human",
  tissue = "global",
  sample_cutoff_method = "reference_guided",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

det$sample_cutoff_summary
```

## PBMC3K Public Smoke Test

`load_pbmc3k_online()` optionally loads the public 10x Genomics PBMC3K dataset.
PBMC3K is useful for small public-data smoke testing. Because PBMCs often have
relatively low mitochondrial fractions and dataset-specific QC structure,
PBMC3K cutoff results should be interpreted as a functionality check rather
than evidence of a universal mitochondrial threshold.

```r
pbmc <- load_pbmc3k_online(method = "10x")
det_pbmc <- SCdetMito(
  pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "human",
  tissue = "PBMC",
  sample_cutoff_method = "reference_guided",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)
```

## Mouse Brain 5k Dataset-Scouting Candidate

`load_mousebrain5k_online()` optionally loads a public 10x Genomics mouse brain
nuclei dataset. It is a real tissue/nuclei candidate for data scouting before
selecting SoftwareX figures. Nuclei datasets may have low mitochondrial
fractions, so visual suitability should be checked.

```r
brain <- load_mousebrain5k_online()
det_brain <- SCdetMito(
  brain,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "mouse",
  tissue = "brain",
  sample_cutoff_method = "reference_guided",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

det_brain$sample_cutoff_summary
```

## Heart10k Solid-Tissue Stress Test

`load_heart10k_online()` optionally loads the public 10x Genomics mouse heart
10k dataset. It is a larger solid-tissue stress-test candidate. It can have
high mitochondrial burden and may produce stringent largest-drop QC decisions,
so reference-guided and first-boundary outputs should be inspected before using
the dataset for SoftwareX figures.

```r
heart <- load_heart10k_online()
det_heart <- SCdetMito(
  heart,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "mouse",
  tissue = "heart",
  sample_cutoff_method = "reference_guided",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

det_heart$sample_cutoff_summary
```

## IFNB Multi-Condition Example

`load_ifnb_online()` optionally loads the public IFNB-stimulated/control PBMC
dataset from SeuratData. It is intended for demonstrating consensus, strictest,
and groupwise multi-sample/group QC strategies. If the dataset is not installed
locally, use `load_ifnb_online(install_if_missing = TRUE)`.

```r
ifnb <- load_ifnb_online()

qc_groupwise <- SCQCmulti(
  ifnb,
  sample_col = "sample",
  group_col = "condition",
  mito_col = "mitoRatio",
  cutoff_strategy = "groupwise",
  scdet_options = list(sample_cutoff_method = "reference_guided"),
  species = "human",
  tissue = "PBMC",
  min_genes = 0,
  min_counts = 0,
  removeDouble = FALSE,
  write_plots = FALSE,
  write_tables = FALSE
)
```

## Dataset Scouting Script

An optional script is provided for comparing public candidates before selecting
SoftwareX figures:

```r
system.file(
  "examples",
  "softwarex_dataset_scouting_online.R",
  package = "SCdetMito"
)
```

The script is not run during package checks. It writes retained-cell and
interval-loss outputs for the built-in demo and any optional public datasets
that can be loaded in the user's environment.

## Main Functions

| Function | Purpose |
| --- | --- |
| `load_demo_seurat()` | Load a lightweight built-in demo object for smoke tests |
| `load_pbmc3k_demo()` | Load a small bundled PBMC3K-derived object for package examples |
| `load_pbmc3k_online()` | Optionally load PBMC3K for small public-data smoke testing |
| `load_mousebrain5k_online()` | Optionally load public 10x mouse brain nuclei data for dataset scouting |
| `load_heart10k_online()` | Optionally load public 10x mouse heart 10k data as a solid-tissue stress test |
| `load_ifnb_online()` | Optionally load IFNB-stimulated/control PBMC data for multi-condition QC examples |
| `SCdetMito_reference_cutoffs()` | Return literature-informed reference cutoffs for interpretation |
| `add_mitoRatio()` | Add mitochondrial fraction to Seurat metadata |
| `SCdetMito()` | Estimate reference-aware adaptive mitochondrial QC cutoffs |
| `SCQCone()` | Apply QC to one Seurat object |
| `SCQCmulti()` | Apply sample-aware multi-sample QC |
| `SCQCbenchmark()` | Compare fixed and adaptive QC strategies |
| `SCdetMito_sensitivity()` | Evaluate cutoff stability across parameter settings |
| `SCdetMito_methods()` | List available detection and cutoff-selection methods |

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
user-adjustable scoring profiles. Default adaptive strategies include
`SCdetMito_largest_drop`, `SCdetMito_first_significant_high`, and
`SCdetMito_reference_guided`. The returned recommendation is a strategy ranked
under the selected scoring profile, not a universal QC rule.

```r
benchmark <- SCQCbenchmark(
  seurat_obj = seu,
  sample_col = "sample",
  group_col = "group",
  species = "human",
  tissue = "global",
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
  sample_cutoff_method = c("largest_drop", "reference_guided"),
  species = "human",
  tissue = "global",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

sensitivity$sensitivity_summary
```

## Interpretation Notes

- SCdetMito provides mitochondrial QC decision support.
- It does not define a universal biological mitochondrial cutoff.
- Reference cutoffs are literature-informed priors, not hard filtering rules.
- Tissue-note reference rows may use species-level priors unless explicitly
  marked as tissue-specific.
- High selected cutoffs should trigger inspection of sample quality, tissue
  handling, post-mortem interval, dissociation stress, disease state, and
  cell-type composition.
- Very stringent largest-drop cutoffs may over-filter real tissue data and
  should be inspected alongside first significant boundaries and reference
  cutoffs.
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
reference cutoff source, sample cutoff summary, and whether cutoffs were
statistically detected or fallback-derived.

## Contact

Silu Hu  
husilu0902@gmail.com  
Repository: https://github.com/ericHuu/SCdetMito
