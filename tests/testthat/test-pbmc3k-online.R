test_that("online real-data loaders are exported", {
  expect_true(is.function(load_pbmc3k_online))
  expect_true(is.function(load_heart10k_online))
  expect_true(is.function(load_ifnb_online))
  expect_true("cache_dir" %in% names(formals(load_pbmc3k_online)))
  expect_true("cache_dir" %in% names(formals(load_heart10k_online)))
  expect_true("group_col" %in% names(formals(load_ifnb_online)))
})

test_that("load_pbmc3k_online can load public data when online tests are enabled", {
  skip_on_cran()
  skip_if_not(Sys.getenv("SCDETMITO_RUN_ONLINE_TESTS") == "true")

  pbmc <- load_pbmc3k_online(method = "10x")
  expect_s4_class(pbmc, "Seurat")
  expect_true("mitoRatio" %in% colnames(pbmc@meta.data))
})

test_that("load_heart10k_online can load public data when online tests are enabled", {
  skip_on_cran()
  skip_if_not(Sys.getenv("SCDETMITO_RUN_ONLINE_TESTS") == "true")

  heart <- load_heart10k_online()
  expect_s4_class(heart, "Seurat")
  expect_true("mitoRatio" %in% colnames(heart@meta.data))
})

test_that("load_ifnb_online can load multi-condition data when online tests are enabled", {
  skip_on_cran()
  skip_if_not(Sys.getenv("SCDETMITO_RUN_ONLINE_TESTS") == "true")

  ifnb <- load_ifnb_online()
  expect_s4_class(ifnb, "Seurat")
  expect_true("condition" %in% colnames(ifnb@meta.data))
  expect_true("sample" %in% colnames(ifnb@meta.data))
  expect_true("mitoRatio" %in% colnames(ifnb@meta.data))
})
