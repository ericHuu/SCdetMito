library(SCdetMito)

# This example downloads the public PBMC3K dataset and therefore requires
# internet access. It is intended as a small public-data smoke test and is not
# run during package checks.

out_dir <- "SCdetMito_pbmc3k_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pbmc <- load_pbmc3k_online(method = "10x")

det <- SCdetMito(
  pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  return_details = TRUE,
  write_tables = TRUE,
  write_plots = TRUE,
  output_dir = out_dir
)

bench <- SCQCbenchmark(
  pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  min_genes = 0,
  min_counts = 0,
  removeDouble = FALSE,
  run_downstream = FALSE,
  write_plots = TRUE,
  output_dir = file.path(out_dir, "benchmark")
)

sens <- SCdetMito_sensitivity(
  pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  return_details = TRUE
)

write.csv(
  sens$sensitivity_summary,
  file.path(out_dir, "pbmc3k_sensitivity_summary.csv"),
  row.names = FALSE
)

sink(file.path(out_dir, "sessionInfo_pbmc3k.txt"))
sessionInfo()
sink()
