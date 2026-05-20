test_that("atlas_benchmark_manifest returns core benchmark candidates", {
  manifest <- atlas_benchmark_manifest()

  expect_true(is.data.frame(manifest))
  expect_true(all(c("dataset_id", "species", "tissue", "role", "journal", "status") %in% colnames(manifest)))
  expect_true("availability_validated" %in% colnames(manifest))
  expect_false(any(manifest$availability_validated))
  expect_true("GSE135893" %in% manifest$dataset_id)
  expect_true("GSE159929" %in% manifest$dataset_id)
  expect_true(any(manifest$species == "chicken"))
  expect_true(any(manifest$species == "duck"))
})

test_that("atlas_benchmark_manifest supports filters", {
  lung_manifest <- atlas_benchmark_manifest(tissue = "lung")
  mouse_manifest <- atlas_benchmark_manifest(species = "mouse")

  expect_true(all(lung_manifest$tissue == "lung"))
  expect_true(all(mouse_manifest$species == "mouse"))
})
