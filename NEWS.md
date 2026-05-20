# SCdetMito 1.3.6

Date: 2026-05-20

SoftwareX submission candidate.

- Added optional online PBMC3K real-data loader through
  `load_pbmc3k_online()`.
- Added optional public 10x mouse heart 10k online loader through
  `load_heart10k_online()` for larger real-tissue demonstration.
- Added optional IFNB-stimulated/control PBMC loader through
  `load_ifnb_online()` for multi-condition QC strategy examples.
- Added optional SoftwareX example scripts for PBMC3K, heart10k, and IFNB
  workflows without bundling large data.
- Clarified that the built-in demo is for smoke tests, PBMC3K is a small
  public-data check, heart10k is the recommended larger real-tissue example,
  and IFNB is the recommended multi-condition workflow example.
- Updated SoftwareX-oriented README, citation metadata, release metadata, and
  author contact information.
- Kept online datasets out of the package to maintain lightweight
  installation.
- Kept standard tests independent of internet access.
- Preserved the core retention-loss enrichment workflow and backward
  compatibility.
- `devtools::test()`: 0 failures, 0 warnings, 6 skipped optional/local tests,
  162 passed expectations.
- `devtools::check(document = FALSE, args = c("--no-manual",
  "--no-build-vignettes"))`: 0 errors, 0 warnings, 1 environment-related NOTE
  for future file timestamp verification.

# SCdetMito 1.3.5

Date: 2026-05-18

- Reframed the core method as retention-loss enrichment-based adaptive
  mitochondrial cutoff detection rather than classical change-point detection.
- Changed the default loss test to `mad_zscore` and the default
  sample-level cutoff selector to `largest_drop`.
- Added explicit support for `p_adj_method = "none"` and a warning for
  adjusted `empirical_tail` sensitivity analyses.
- Aligned `SCQCmulti(cutoff_strategy = "consensus")` with the documented
  all-sample sample-supported global cutoff, even when `group_by` is supplied.
- Added explicit retained-cell profile, interval-loss, candidate-filter,
  adjusted p-value, significant-interval, and sample-cutoff summary outputs.
- Added `SCdetMito_sensitivity()` for parameter robustness analyses across
  bin widths, alpha levels, loss tests, and sample cutoff selection rules.
- Added sample-supported global cutoff metadata and clearer consensus,
  strictest, and groupwise aggregation outputs.
- Added QC provenance tracking under `seurat_obj@misc$SCdetMito_QC` and safer
  optional DoubletFinder handling.
- Updated benchmark outputs with fixed 5/10/15/20 percent baselines, adaptive
  SCdetMito strategies, built-in scoring profiles, and manuscript-quality
  plotting filenames.
- Hardened `add_mitoRatio()` to stop clearly when no mitochondrial features are
  detected and to record mitochondrial feature metadata in `@misc$SCdetMito`.

# SCdetMito 1.3.3

Date: 2026-05-14

- Added package-installed documentation under `inst/doc/`, including a compact
  HTML user guide for offline package use.
- Added package citation metadata in `inst/CITATION`.
- Tightened package metadata for a Seurat-oriented workflow by requiring
  `R >= 4.1.0`.
- Clarified publication readiness, package-versus-analysis file layout, and
  manuscript evidence gaps.
- Kept public benchmark datasets and manuscript-scale analyses outside the
  installable R package while retaining the lightweight PBMC3K demo object.

# SCdetMito 1.3.2

Date: 2026-05-14

- Moved `DoubletFinder` from a hard dependency to an optional suggested
  dependency and changed `removeDouble` defaults to `FALSE`.
- Added the robust `mad_zscore` loss-enrichment test option.
- Removed the local absolute-path installation call from `tests/testthat.R`.
- Added a strict reviewer-style Chinese method assessment to guide manuscript
  and package hardening.
- Prepared the package for lightweight PBMC3K-based examples and source
  package release.

# SCdetMito 1.3.1

Date: 2026-05-08

- Added centralized QC threshold validation for gene, count, and mitochondrial
  bounds.
- Allowed numeric mitochondrial cutoff inputs from `1` to `100` to be treated
  as percentages and converted to fractions with an explicit warning.
- Added tests for percentage-style mitochondrial cutoff inputs and invalid
  mitochondrial bound combinations.
- Standardized R source, test, and vignette formatting.
- Refreshed the local HTML documentation workflow.

# SCdetMito 1.3.0

Date: 2026-04-19

- Added group-aware multi-sample mitochondrial QC strategies.
- Added benchmark helpers for comparing fixed and adaptive QC strategies.
- Added bundled multi-group Seurat demo data.
- Added cross-species mitochondrial feature presets for livestock and birds.
