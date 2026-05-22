test_that("legacy PBMC3K demo loader remains a compatibility wrapper", {
  expect_warning(
    seu <- load_pbmc3k_demo(),
    "legacy alias"
  )

  expect_s4_class(seu, "Seurat")
  expect_equal(ncol(seu), ncol(load_demo_pbmc()))
  expect_true(all(c("sample", "group", "mitoRatio") %in% colnames(seu@meta.data)))
})
