# SCdetMito 1.4.3

Date: 2026-05-23

## Provenance clarification and documentation harmonization

- Fixed `SCQCmulti()` groupwise provenance so top-level
  `recommended_cutoff` and `cutoff_applied` are consistent with the final
  applied group-level cutoff.
- Added explicit `sample_cutoff_summary` and `group_cutoff_summary` aliases in
  `SCQCmulti()` provenance.
- Clarified the distinction among sample-level cutoffs, group-level cutoffs,
  and final applied cutoffs.
- Harmonized wording across README, vignettes, roxygen documentation, and
  examples.
- Standardized scientific terminology to describe SCdetMito as a
  reference-aware and sample-aware mitochondrial QC decision-support framework.
- Avoided language implying universal, true, or biologically optimal
  mitochondrial thresholds.

# SCdetMito 1.4.2

Date: 2026-05-22

## Release candidate bug fixes and documentation cleanup

- Removed IFNB from core package examples and exports because processed IFNB
  objects may lack mitochondrial genes or a valid mitochondrial fraction
  column.
- Added explicit `recommended_cutoff` reporting independent of the
  user-selected cutoff mode.
- Added recommendation metadata fields including `recommended_method`,
  `recommended_reason`, `recommendation_level`, `recommendation_warning`, and
  `recommendation_source`.
- Clarified the relationship among `selected_cutoff`, `detected_cutoff`, and
  `recommended_cutoff`.
- Updated `SCQCone()` and `SCQCmulti()` to apply `recommended_cutoff` by
  default when adaptive filtering is requested, with provenance recording the
  applied cutoff source.
- Added recommendation and reference summary CSV outputs when
  `write_tables = TRUE`.
- Internal demo strategy simplified to one practical multi-donor PBMC demo.
- Added bundled demo data derived from public 10x Genomics GEM-X Human PBMC
  Donor 1-4 datasets.
- Updated `load_demo_pbmc()` as the primary internal demo loader and retained
  `load_demo_seurat()` as a backward-compatible alias.
- Updated examples, README, and vignettes to use the internet-independent
  multi-donor PBMC demo for sample-aware and group-aware workflows.
- Clarified that muscle and other tissue datasets are handled in separate
  manuscript-level validation workflows outside the core package.
- Fixed mitoRatio consistency between precomputed and auto-computed workflows.
- Clarified mitochondrial ratio overwrite/recompute policy.
- Consolidated mitochondrial ratio handling through `ensure_mito_ratio()`.
- Refactored `add_mitoRatio()` as a backward-compatible wrapper.
- Removed or consolidated duplicated, unused, or conflicting helper functions.
- Expanded mitochondrial exact feature fallback presets.
- Updated README, vignettes, examples, and pkgdown documentation for SoftwareX
  readiness.

# SCdetMito 1.4.0

Date: 2026-05-20

## Reference-aware mitochondrial QC decision support

- Added literature-informed reference cutoff table through
  `SCdetMito_reference_cutoffs()`.
- Added reference-aware cutoff reporting, including reference cutoff, first
  significant high-to-low boundary, largest-drop cutoff, and selected cutoff.
- Added new cutoff selection modes including `first_significant_high`,
  `first_significant_low`, `reference_guided`, and
  `largest_drop_with_reference_guard`.
- Added reference-deviation, low-retention, and fallback warning fields in
  `sample_cutoff_summary`.
- Updated plots to distinguish reference, first significant, largest-drop, and
  selected cutoffs where available.
- Updated provenance tracking for reference-aware QC decisions.
- Updated benchmark defaults to include largest-drop, first-boundary, and
  reference-guided SCdetMito adaptive strategies.
- Added `ensure_mito_ratio()` as the unified mitochondrial fraction preparation
  function.
- Refactored `add_mitoRatio()` as a backward-compatible wrapper around
  `ensure_mito_ratio()`.
- Updated `SCdetMito()`, `SCQCone()`, `SCQCmulti()`, and `SCQCbenchmark()` to
  share the same mitochondrial ratio validation and calculation workflow.
- Updated optional public data loaders to use the shared mitochondrial ratio
  preparation workflow.
- Consolidated duplicated mitochondrial feature detection logic.
- Removed incomplete manual mitochondrial gene lists from documentation
  examples.
- Updated README and documentation to clarify that SCdetMito provides QC
  decision support rather than universal biological mitochondrial thresholds.
- Corrected reference metadata and DOI information.
- Preserved backward compatibility with `detected_cutoff`, legacy output files,
  and legacy result fields where possible.

# SCdetMito 1.3.6

Date: 2026-05-20

SoftwareX submission candidate.

- Added optional public PBMC3K online loader for small real-data smoke testing.
- Added optional public 10x mouse heart 10k online loader as a larger
  solid-tissue stress-test candidate.
- Added optional public mouse brain 5k nuclei loader as a real tissue/nuclei
  candidate for dataset scouting.
- Improved optional online loader handling for public-data examples.
- Added optional dataset-scouting example script to compare public datasets
  before selecting SoftwareX figures.
- Clarified that PBMC3K is a small public-data check, heart10k is a
  solid-tissue stress-test candidate, and mouse brain 5k is a tissue/nuclei
  scouting candidate.
- Updated SoftwareX-oriented README, citation metadata, release metadata, and
  author contact information.
- Kept online datasets out of the package to maintain lightweight
  installation.
- Kept standard tests independent of internet access.
- Preserved the core retention-loss enrichment workflow and backward
  compatibility.

# SCdetMito 1.3.5

Date: 2026-05-18

- Reframed the core method as retention-loss enrichment-based adaptive
  mitochondrial cutoff detection rather than a breakpoint model.
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
