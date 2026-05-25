# SCdetMito

SCdetMito is an R package for reference-aware and sample-aware adaptive
mitochondrial quality control in single-cell RNA-seq data.

Version: `1.4.3`  
Release date: `2026-05-23`

## Overview

SCdetMito provides mitochondrial QC decision support for Seurat workflows. It
constructs retained-cell profiles across candidate mitochondrial cutoffs,
quantifies interval-specific cell loss, reports multiple interpretable cutoff
candidates, integrates literature-informed reference cutoffs, and supports
multi-sample QC strategies, benchmarking, sensitivity analysis, cutoff
confidence reporting, warning fields, and provenance tracking.

The package does not define a universal biological mitochondrial cutoff.
Detected cutoffs should be interpreted with retained-cell profiles,
interval-loss profiles, reference values, cutoff source, cutoff confidence, and
study context.

SCdetMito keeps internal demo data limited to a practical multi-sample PBMC
example. Larger multi-tissue and multi-species validation workflows are
maintained outside the core package.

## Installation

```r
install.packages("remotes")
remotes::install_github("ericHuu/SCdetMito@v1.4.3", build_vignettes = FALSE)
```

For local source installation:

```r
install.packages(
  "/path/to/SCdetMito_1.4.3.tar.gz",
  repos = NULL,
  type = "source"
)
```

## Quick Start

The bundled demo is a downsampled multi-donor PBMC Seurat object derived from
public 10x Genomics GEM-X Human PBMC Donor 1-4 datasets. It contains four
donor-level samples, two demographic groups, raw counts, retained mitochondrial
features, and a precomputed `mitoRatio` column calculated with
`ensure_mito_ratio()`.

```r
library(SCdetMito)

pbmc <- load_demo_pbmc()

det <- SCdetMito(
  pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "human",
  tissue = "PBMC",
  sample_cutoff_method = "reference_guided",
  return_details = TRUE
)

det$sample_cutoff_summary
```

`load_demo_seurat()` is retained as a backward-compatible alias for
`load_demo_pbmc()`.

## Multi-Sample Example

```r
qc_groupwise <- SCQCmulti(
  pbmc,
  sample_col = "sample",
  group_col = "group",
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

The built-in demo is intended for examples, software tests, and sample-aware or
group-aware workflow checks. It is not intended for biological inference.
When adaptive mitochondrial filtering is requested, `SCQCone()` and
`SCQCmulti()` apply `recommended_cutoff` by default and record the applied
cutoff source in object provenance. Set `use_recommended_cutoff = FALSE` to
apply the user-selected `selected_cutoff` instead.

`SCQCmulti()` stores both sample-level cutoff evidence and strategy-level
applied cutoffs. The field `cutoff_applied` records the actual filtering
cutoff. `sample_cutoff_summary` retains per-sample evidence, and
`group_cutoff_summary` records group-level aggregation when `cutoff_strategy =
"groupwise"`.

## Reference-Aware Cutoff Interpretation

SCdetMito separates several cutoff concepts that can differ in real datasets:

| Cutoff type | Meaning |
| --- | --- |
| `reference_cutoff` | Literature-informed species/tissue prior for interpretation |
| `first_significant_cutoff_high` | Most permissive significant retention-loss boundary when scanning from high to low cutoff |
| `largest_drop_cutoff` | Significant cutoff with the largest interval-specific cell loss |
| `selected_cutoff` | Cutoff selected under `sample_cutoff_method` |
| `detected_cutoff` | Backward-compatible alias of `selected_cutoff` |
| `recommended_cutoff` | SCdetMito recommendation from the default reference-aware policy |
| `cutoff_applied` | Actual cutoff used for filtering in `SCQCone()` or `SCQCmulti()` |
| `sample_cutoff_summary` | Per-sample cutoff evidence and recommendation table |
| `group_cutoff_summary` | Group-level aggregation table for groupwise multi-sample QC |

Reference cutoffs are not hard filtering rules. They help flag cases where a
data-driven cutoff is far above common literature expectations or where a
largest-drop cutoff is very stringent and retains few cells. When no exact
tissue-specific value is available, SCdetMito falls back to a species-level
reference and reports this explicitly.

Even when users select a strict or exploratory cutoff mode, SCdetMito reports a
`recommended_cutoff` using the reference-aware policy. The user-selected cutoff
is preserved as `selected_cutoff`, and `detected_cutoff` remains a
backward-compatible alias of `selected_cutoff`.

```r
det <- SCdetMito(
  pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "human",
  tissue = "PBMC",
  sample_cutoff_method = "largest_drop",
  return_details = TRUE
)

det$sample_cutoff_summary[, c(
  "sample",
  "selected_cutoff",
  "recommended_cutoff",
  "reference_cutoff",
  "first_significant_cutoff_high",
  "largest_drop_cutoff",
  "recommendation_level",
  "recommended_reason"
)]
```

## Mitochondrial Ratio Handling

SCdetMito expects mitochondrial content as a fraction between 0 and 1. If an
existing mitochondrial column is present, `ensure_mito_ratio()` validates that
it is numeric and converts percent-scale values to fractions when possible. If
the column is absent, SCdetMito can calculate it from retained mitochondrial
features using species-aware gene-symbol patterns and feature metadata when
available.

Existing mitochondrial columns are preserved by default after validation.
Recalculation occurs only when the column is absent or when
`recompute_mito = TRUE` is supplied.

```r
pbmc <- ensure_mito_ratio(
  pbmc,
  mito_col = "mitoRatio",
  species = "human"
)
```

Users can provide a complete mitochondrial feature vector when automatic
detection is not appropriate:

```r
mito_features <- grep("^MT-", rownames(pbmc), value = TRUE)

det <- SCdetMito(
  pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  mito_features = mito_features,
  auto_add_mito = TRUE,
  recompute_mito = FALSE
)
```

Manual mitochondrial feature lists should include all mitochondrial features
retained in the object. Short illustrative subsets should not be used for real
QC. If mitochondrial genes were removed during preprocessing, mitochondrial
fractions cannot be reconstructed from the Seurat object alone; provide a
precomputed `mito_col` or raw counts retaining mitochondrial genes.

## Main Functions

| Function | Purpose |
| --- | --- |
| `load_demo_pbmc()` | Load the bundled multi-donor public PBMC demo object |
| `load_demo_seurat()` | Backward-compatible alias for `load_demo_pbmc()` |
| `SCdetMito()` | Estimate reference-aware adaptive mitochondrial QC cutoffs |
| `ensure_mito_ratio()` | Validate or calculate mitochondrial fraction metadata |
| `add_mitoRatio()` | Backward-compatible wrapper for mitochondrial fraction calculation |
| `SCdetMito_reference_cutoffs()` | Return literature-informed reference cutoffs for interpretation |
| `SCQCone()` | Apply QC to one Seurat object |
| `SCQCmulti()` | Apply sample-aware multi-sample QC |
| `SCQCbenchmark()` | Compare fixed and adaptive QC strategies |
| `SCdetMito_sensitivity()` | Evaluate cutoff stability across parameter settings |
| `SCdetMito_methods()` | List available detection and cutoff-selection methods |

Optional online public-data loaders remain available as user-run utilities, but
they are not used by package tests, vignettes, or the core workflow examples.

Advanced and backward-compatible utilities remain exported for scripted
workflows:

| Function | Purpose |
| --- | --- |
| `check_seu()` | Validate Seurat metadata columns and optionally normalize fraction-scale metadata |
| `SCQC_processedPlots()` | Standalone QC plotting helper used by legacy and scripted workflows |
| `atlas_benchmark_manifest()` | Reference-only planning manifest; it does not validate dataset availability |

## Benchmarking and Sensitivity

```r
benchmark <- SCQCbenchmark(
  seurat_obj = pbmc,
  sample_col = "sample",
  group_col = "group",
  species = "human",
  tissue = "PBMC",
  min_genes = 0,
  min_counts = 0,
  removeDouble = FALSE,
  run_downstream = FALSE,
  plot = FALSE,
  output_dir = tempdir()
)

benchmark$recommended_strategies
```

`SCQCbenchmark()` returns strategies ranked under the selected scoring profile.
It does not identify a universally applicable QC strategy.

```r
sensitivity <- SCdetMito_sensitivity(
  seurat_obj = pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  sample_cutoff_method = c("largest_drop", "reference_guided"),
  species = "human",
  tissue = "PBMC",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

sensitivity$sensitivity_summary
```

## Optional Public Data Workflows

SCdetMito is not intended to be a public dataset loader collection. Optional
online loaders are retained as convenience utilities for PBMC3K, mouse heart
10k, and mouse brain 5k public examples, but they are not part of the core
package demo workflow. Core examples use the bundled multi-donor PBMC demo and
do not require internet access.

Processed public objects may not retain mitochondrial genes or a valid
mitochondrial fraction column. Run SCdetMito on such objects only after
providing raw counts with mitochondrial features or a valid precomputed
`mito_col`.

Manuscript-scale public data analyses should be maintained outside the core
package. Optional public-data scripts are installed under `inst/examples/`,
including a dataset-scouting script:

```r
system.file(
  "examples",
  "softwarex_dataset_scouting_online.R",
  package = "SCdetMito"
)
```

## Interpretation Notes

- SCdetMito provides mitochondrial QC decision support.
- It does not define a universal biological mitochondrial cutoff.
- Reference cutoffs are literature-informed priors, not hard filtering rules.
- Very stringent largest-drop cutoffs may over-filter some tissue datasets and
  should be inspected alongside first significant boundaries and reference
  cutoffs.
- Fallback-derived cutoffs are explicitly labeled and should be interpreted
  cautiously.
- Mitochondrial filtering should be combined with other QC metrics such as
  `nFeature`, `nCount`, doublet detection, ambient RNA assessment, cell type
  composition, and biological context.
- Benchmark scores are user-adjustable summaries and do not imply a
  universally applicable QC strategy.

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
