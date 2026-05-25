library(SCdetMito)

# This example downloads the public 10x Genomics mouse heart 10k dataset.
# It requires internet access and is not run during package checks.
# The dataset is intended as a larger solid-tissue stress-test candidate for
# SoftwareX data scouting.

out_dir <- "SCdetMito_heart10k_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

heart <- load_heart10k_online()

det <- SCdetMito(
  heart,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "mouse",
  tissue = "heart",
  sample_cutoff_method = "reference_guided",
  return_details = TRUE,
  write_tables = TRUE,
  write_plots = TRUE,
  output_dir = out_dir
)

bench <- SCQCbenchmark(
  heart,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "mouse",
  tissue = "heart",
  min_genes = 0,
  min_counts = 0,
  removeDouble = FALSE,
  run_downstream = FALSE,
  write_plots = TRUE,
  output_dir = file.path(out_dir, "benchmark")
)

sens <- SCdetMito_sensitivity(
  heart,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "mouse",
  tissue = "heart",
  return_details = TRUE
)

write.csv(
  sens$sensitivity_summary,
  file.path(out_dir, "heart10k_sensitivity_summary.csv"),
  row.names = FALSE
)

sink(file.path(out_dir, "sessionInfo_heart10k.txt"))
sessionInfo()
sink()
