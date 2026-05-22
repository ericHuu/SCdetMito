library(SCdetMito)

# This optional script compares candidate public datasets before selecting
# SoftwareX figures. It downloads or loads public data only when run manually
# and is not used during package checks.

out_dir <- "SCdetMito_dataset_scouting"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

run_one <- function(name, obj, species = NULL, tissue = NULL) {
  message("Running SCdetMito on: ", name)
  ds_dir <- file.path(out_dir, name)
  dir.create(ds_dir, showWarnings = FALSE, recursive = TRUE)

  det <- SCdetMito(
    obj,
    sample_col = "sample",
    mito_col = "mitoRatio",
    species = species,
    tissue = tissue,
    sample_cutoff_method = "reference_guided",
    return_details = TRUE,
    write_tables = TRUE,
    write_plots = TRUE,
    output_dir = ds_dir
  )

  summary <- det$sample_cutoff_summary
  summary$dataset_name <- name
  summary
}

results <- list()

seu <- load_demo_seurat()
results$demo <- run_one("demo", seu, species = "human", tissue = "global")

try({
  pbmc <- load_pbmc3k_online(method = "10x")
  results$pbmc3k <- run_one("pbmc3k", pbmc, species = "human", tissue = "PBMC")
}, silent = TRUE)

try({
  heart <- load_heart10k_online()
  results$heart10k <- run_one("heart10k", heart, species = "mouse", tissue = "heart")
}, silent = TRUE)

try({
  brain <- load_mousebrain5k_online()
  results$mousebrain5k <- run_one("mousebrain5k", brain, species = "mouse", tissue = "brain")
}, silent = TRUE)

scout_summary <- do.call(rbind, results)

write.csv(
  scout_summary,
  file.path(out_dir, "SCdetMito_dataset_scouting_summary.csv"),
  row.names = FALSE
)

print(scout_summary)

# Suggested suitability review:
# - fallback_used is FALSE
# - cutoff_confidence is moderate or high
# - retention_fraction_at_cutoff is not extremely low
# - selected_cutoff is biologically interpretable
# - reference and warning fields support transparent interpretation
# - interval-loss plots are visually clear

sink(file.path(out_dir, "sessionInfo_dataset_scouting.txt"))
sessionInfo()
sink()
