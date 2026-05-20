test_that("lung IPF GEO object can run group-aware workflows when present locally", {
  # Optional local validation test; package checks do not require this file.
  project_root <- normalizePath(
    file.path(testthat::test_path("..", "..")),
    winslash = "/",
    mustWork = TRUE
  )
  public_path <- file.path(project_root, "data-public", "lung_ipf_geo", "GSE157376_lung_ipf_seurat.rds")
  testthat::skip_if_not(file.exists(public_path))

  seu <- readRDS(public_path)
  output_dir <- tempfile("scdetmito-lung-ipf-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  qc_groupwise <- suppressWarnings(SCQCmulti(
    seu,
    by = "sample",
    group_by = "group",
    mode = "all",
    cutoff_strategy = "groupwise",
    max_mito = "SCdetMito",
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = TRUE,
    output_dir = output_dir
  ))

  expect_s4_class(qc_groupwise, "Seurat")
  expect_true("group" %in% colnames(qc_groupwise@meta.data))
  expect_true(file.exists(file.path(output_dir, "SCQCmulti_group_cutoffs.csv")))
})
