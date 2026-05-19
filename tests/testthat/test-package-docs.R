test_that("installed package documentation and citation files are present", {
  project_file <- function(...) {
    file.path(testthat::test_path("..", ".."), ...)
  }

  package_file <- function(...) {
    installed_path <- system.file(..., package = "SCdetMito")
    if (nzchar(installed_path) && file.exists(installed_path)) {
      return(installed_path)
    }

    file.path(testthat::test_path("..", ".."), "inst", ...)
  }

  expect_true(
    file.exists(package_file("doc", "SCdetMito_package_docs.html")) ||
      file.exists(project_file("vignettes", "SCdetMito_package_docs.Rmd"))
  )
  expect_true(
    file.exists(package_file("doc", "SCdetMito_package_docs.Rmd")) ||
      file.exists(package_file("doc", "SCdetMito_package_docs.md")) ||
      file.exists(project_file("vignettes", "SCdetMito_package_docs.Rmd"))
  )
  expect_true(file.exists(package_file("CITATION")))
})
