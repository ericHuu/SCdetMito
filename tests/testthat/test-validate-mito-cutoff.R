test_that("validate_mito_cutoff reports supportive secondary evidence", {
  counts <- matrix(
    rep(c(5, 3, 2, 1), 30),
    nrow = 4,
    dimnames = list(paste0("g", 1:4), paste0("c", 1:30))
  )
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts))
  seu$sample <- rep(c("A", "B"), each = 15)
  seu$mitoRatio <- c(rep(0.05, 20), rep(0.20, 10))
  seu$nFeature_RNA <- c(rep(1000, 20), rep(400, 10))
  seu$nCount_RNA <- c(rep(3000, 20), rep(1000, 10))

  result <- validate_mito_cutoff(
    seu,
    cutoff = 0.10,
    sample_col = "sample",
    min_cells_per_group = 2
  )

  expect_s3_class(result, "data.frame")
  expect_equal(result$evidence_status[result$scope == "overall"], "supportive")
  expect_true(all(c(
    "fraction_at_or_below",
    "nFeature_relative_change",
    "nCount_relative_change",
    "complexity_median_above"
  ) %in% colnames(result)))
  expect_equal(nrow(result), 3)
})

test_that("validate_mito_cutoff never chooses or changes the supplied cutoff", {
  counts <- matrix(
    rep(c(5, 3, 2, 1), 10),
    nrow = 4,
    dimnames = list(paste0("g", 1:4), paste0("c", 1:10))
  )
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts))
  seu$mitoRatio <- seq(0.01, 0.20, length.out = 10)

  result <- validate_mito_cutoff(
    seu,
    cutoff = 10,
    cutoff_scale = "percent",
    min_cells_per_group = 2
  )

  expect_equal(result$cutoff, 0.10)
  expect_equal(result$evidence_status, "inconclusive")
  expect_false(any(grepl("recommended|selected", colnames(result))))
})

test_that("validate_mito_cutoff validates required metadata", {
  counts <- matrix(
    1,
    nrow = 3,
    ncol = 6,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:6))
  )
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts))

  expect_error(
    validate_mito_cutoff(seu, cutoff = 0.10, mito_col = "missing"),
    "Missing required metadata"
  )
})
