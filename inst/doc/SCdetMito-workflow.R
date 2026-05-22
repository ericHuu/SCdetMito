## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(SCdetMito)

pbmc <- load_demo_pbmc()
table(pbmc$sample)
table(pbmc$group)
summary(pbmc$mitoRatio)

## -----------------------------------------------------------------------------
det <- SCdetMito(
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

det$sample_cutoff_summary[, c(
  "sample",
  "reference_cutoff",
  "first_significant_cutoff_high",
  "largest_drop_cutoff",
  "selected_cutoff",
  "detected_cutoff",
  "warning_level"
)]

## -----------------------------------------------------------------------------
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

qc_groupwise@misc$SCdetMito_QC$group_plan

## -----------------------------------------------------------------------------
benchmark <- SCQCbenchmark(
  pbmc,
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

## -----------------------------------------------------------------------------
sensitivity <- SCdetMito_sensitivity(
  pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "human",
  tissue = "PBMC",
  sample_cutoff_method = c("largest_drop", "reference_guided"),
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

sensitivity$sensitivity_summary

