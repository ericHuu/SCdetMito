library(SCdetMito)

# This optional example loads the public IFNB-stimulated/control PBMC dataset
# from SeuratData. It requires SeuratData and is not run during package checks.
# Processed IFNB objects may not retain mitochondrial genes or mitoRatio, so
# mitochondrial QC steps are guarded.

out_dir <- "SCdetMito_ifnb_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ifnb <- load_ifnb_online(add_mito = FALSE)

if ("mitoRatio" %in% colnames(ifnb@meta.data)) {
  det <- SCdetMito(
    ifnb,
    sample_col = "sample",
    mito_col = "mitoRatio",
    species = "human",
    tissue = "PBMC",
    sample_cutoff_method = "reference_guided",
    return_details = TRUE,
    write_tables = TRUE,
    write_plots = TRUE,
    output_dir = file.path(out_dir, "SCdetMito")
  )

  qc_groupwise <- SCQCmulti(
    ifnb,
    sample_col = "sample",
    group_col = "condition",
    mito_col = "mitoRatio",
    cutoff_strategy = "groupwise",
    species = "human",
    tissue = "PBMC",
    scdet_options = list(sample_cutoff_method = "reference_guided"),
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    write_plots = FALSE,
    write_tables = TRUE,
    output_dir = file.path(out_dir, "groupwise")
  )

  bench <- SCQCbenchmark(
    ifnb,
    sample_col = "sample",
    group_col = "condition",
    mito_col = "mitoRatio",
    species = "human",
    tissue = "PBMC",
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    run_downstream = FALSE,
    write_plots = TRUE,
    output_dir = file.path(out_dir, "benchmark")
  )
} else {
  message(
    "Skipping IFNB mitochondrial QC: no mitoRatio column is available. ",
    "Provide raw counts with mitochondrial genes or a valid mito_col before ",
    "running SCdetMito or SCQCmulti."
  )
}

sink(file.path(out_dir, "sessionInfo_ifnb.txt"))
sessionInfo()
sink()
