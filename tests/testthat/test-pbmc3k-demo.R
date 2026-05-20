test_that("bundled lightweight PBMC3K object supports package workflows", {
  seu <- load_pbmc3k_demo()

  expect_s4_class(seu, "Seurat")
  expect_equal(ncol(seu), 300)
  expect_true(all(c("sample", "group", "mitoRatio") %in% colnames(seu@meta.data)))
  expect_true(length(unique(seu@meta.data$sample)) >= 2)

  cutoff <- suppressWarnings(SCdetMito(
    seu,
    by = "sample",
    mitoRatio = "mitoRatio",
    min_cut = 0.01,
    max_cut = 0.25,
    plot = FALSE,
    table_out = FALSE
  ))
  expect_true(is.numeric(cutoff))
  expect_true(cutoff >= 0.01 && cutoff <= 0.25)

  qc_obj <- suppressWarnings(SCQCmulti(
    seu,
    by = "sample",
    mode = "all",
    cutoff_strategy = "consensus",
    max_mito = "SCdetMito",
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = FALSE
  ))

  expect_s4_class(qc_obj, "Seurat")
  expect_true(ncol(qc_obj) <= ncol(seu))
  expect_true(!is.null(qc_obj@misc$SCdetMito_QC$sample_plan))
})
