test_that("public 10x PBMC integration object can run the main workflows", {
  # Optional local validation test; package checks do not require this file.
  project_root <- normalizePath(
    file.path(testthat::test_path("..", "..")),
    winslash = "/",
    mustWork = TRUE
  )
  public_path <- file.path(project_root, "data-public", "pbmc_10x", "pbmc3k_pbmc4k_public_seurat.rds")
  testthat::skip_if_not(file.exists(public_path))

  seu <- readRDS(public_path)
  output_dir <- tempfile("scdetmito-public-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  qc_single <- suppressWarnings(SCQCone(
    subset(seu, subset = sample == "pbmc3k"),
    max_mito = "SCdetMito",
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = FALSE,
    output_dir = output_dir
  ))

  qc_multi <- suppressWarnings(SCQCmulti(
    seu,
    by = "sample",
    mode = "all",
    cutoff_strategy = "consensus",
    max_mito = "SCdetMito",
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = TRUE,
    output_dir = output_dir
  ))

  expect_s4_class(qc_single, "Seurat")
  expect_s4_class(qc_multi, "Seurat")
  expect_true("mitoRatio" %in% colnames(qc_multi@meta.data))
  expect_true(file.exists(file.path(output_dir, "SCQCmulti_cutoff_plan.csv")))
})
