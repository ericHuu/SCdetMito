## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(SCdetMito)

seu <- load_demo_seurat()
seu
table(seu$group)
table(seu$sample)

## -----------------------------------------------------------------------------
sample_s1 <- subset(seu, subset = sample == "S1")

qc_single <- SCQCone(
  seurat_obj = sample_s1,
  min_genes = 0,
  min_counts = 0,
  max_mito = 0.10,
  removeDouble = FALSE,
  plot = FALSE,
  table_out = FALSE
)

ncol(sample_s1)
ncol(qc_single)

## -----------------------------------------------------------------------------
qc_groupwise <- SCQCmulti(
  seurat_obj = seu,
  by = "sample",
  group_by = "group",
  mode = "all",
  cutoff_strategy = "groupwise",
  species = "human",
  tissue = "global",
  scdet_options = list(sample_cutoff_method = "reference_guided"),
  min_genes = 0,
  min_counts = 0,
  max_mito = "SCdetMito",
  removeDouble = FALSE,
  plot = FALSE,
  table_out = FALSE
)

ncol(seu)
ncol(qc_groupwise)
qc_groupwise@misc$SCdetMito_QC$group_plan

## -----------------------------------------------------------------------------
benchmark <- SCQCbenchmark(
  seurat_obj = seu,
  sample_col = "sample",
  group_col = "group",
  strategies = list(
    list(name = "fixed_10", mode = "all", cutoff_strategy = "consensus", max_mito = 0.10),
    list(name = "fixed_15", mode = "all", cutoff_strategy = "consensus", max_mito = 0.15),
    list(
      name = "SCdetMito_reference_guided",
      mode = "all",
      cutoff_strategy = "groupwise",
      max_mito = "SCdetMito",
      scdet_options = list(sample_cutoff_method = "reference_guided")
    )
  ),
  species = "human",
  tissue = "global",
  min_genes = 0,
  min_counts = 0,
  removeDouble = FALSE,
  run_downstream = FALSE,
  plot = FALSE,
  save_objects = FALSE,
  output_dir = tempdir()
)

benchmark$recommendation
benchmark$recommended_strategies
benchmark$scores[, c("strategy", "overall_score", "retention_rate", "median_mito_after")]

