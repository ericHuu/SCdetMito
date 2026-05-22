test_that("online real-data loaders are exported", {
  expect_true(is.function(load_pbmc3k_online))
  expect_true(is.function(load_heart10k_online))
  expect_true(is.function(load_mousebrain5k_online))
  expect_true("cache_dir" %in% names(formals(load_pbmc3k_online)))
  expect_true("cache_dir" %in% names(formals(load_heart10k_online)))
  expect_true("cache_dir" %in% names(formals(load_mousebrain5k_online)))
  expect_true("url" %in% names(formals(load_heart10k_online)))
  expect_true("local_file" %in% names(formals(load_heart10k_online)))
  expect_true("url" %in% names(formals(load_mousebrain5k_online)))
  expect_true("local_file" %in% names(formals(load_mousebrain5k_online)))
  expect_false("load_ifnb_online" %in% getNamespaceExports("SCdetMito"))
})

test_that("README and vignettes do not present IFNB workflows", {
  package_root <- file.path(testthat::test_path("..", ".."))
  doc_files <- c(
    file.path(package_root, "README.md"),
    list.files(
      file.path(package_root, "vignettes"),
      pattern = "[.]Rmd$",
      full.names = TRUE
    )
  )
  testthat::skip_if_not(all(file.exists(doc_files)), "Source README/vignettes are unavailable in this test context.")
  doc_text <- paste(vapply(doc_files, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }, character(1)), collapse = "\n")

  expect_false(grepl("load_ifnb_online|\\bIFNB\\b|\\bifnb\\b", doc_text))
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

test_that("load_mousebrain5k_online can load public data when online tests are enabled", {
  skip_on_cran()
  skip_if_not(Sys.getenv("SCDETMITO_RUN_ONLINE_TESTS") == "true")

  brain <- load_mousebrain5k_online()
  expect_s4_class(brain, "Seurat")
  expect_true("mitoRatio" %in% colnames(brain@meta.data))
})
