# Reproducible offline quick start using the bundled multi-donor PBMC object.
library(SCdetMito)

pbmc <- load_demo_pbmc()

detection <- SCdetMito(
  pbmc,
  sample_col = "sample",
  mito_col = "mitoRatio",
  species = "human",
  tissue = "PBMC",
  return_details = TRUE,
  write_plots = FALSE,
  write_tables = FALSE
)

print(detection$sample_cutoff_summary[, c(
  "sample",
  "reference_cutoff",
  "data_candidate_cutoff",
  "recommended_cutoff",
  "recommendation_status",
  "recommended_auto_apply_eligible"
)])

if (is.finite(detection$recommended_cutoff)) {
  print(validate_mito_cutoff(
    pbmc,
    cutoff = detection$recommended_cutoff,
    mito_col = "mitoRatio",
    sample_col = "sample"
  ))
}
