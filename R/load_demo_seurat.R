#' Load the bundled multi-donor PBMC demo object
#'
#' @description
#' `load_demo_pbmc()` loads the bundled SCdetMito demonstration object derived
#' from public 10x Genomics GEM-X Human PBMC donor datasets. The object contains
#' four donors, sample and group metadata, retained raw counts, mitochondrial
#' features, and a precomputed `mitoRatio` column calculated with
#' [ensure_mito_ratio()]. It is intended for software examples and tests, not
#' biological inference.
#'
#' @return A Seurat object stored in
#'   `inst/extdata/demo_pbmc_multidonor_seurat.rds`.
#' @export
load_demo_pbmc <- function() {
  demo_path <- system.file("extdata", "demo_pbmc_multidonor_seurat.rds", package = "SCdetMito")
  if (!nzchar(demo_path)) {
    stop("The bundled multi-donor PBMC demo Seurat object could not be located.", call. = FALSE)
  }

  readRDS(demo_path)
}

#' Load the bundled Seurat demo object
#'
#' @description
#' `load_demo_seurat()` is retained for backward compatibility and now calls
#' [load_demo_pbmc()]. New examples should use `load_demo_pbmc()` directly.
#'
#' @return A Seurat object from [load_demo_pbmc()].
#' @export
load_demo_seurat <- function() {
  load_demo_pbmc()
}

#' Legacy PBMC3K demo loader
#'
#' @description
#' `load_pbmc3k_demo()` is retained as a backward-compatible legacy alias. The
#' current bundled demo is a multi-donor public GEM-X PBMC object loaded by
#' [load_demo_pbmc()], not the earlier lightweight PBMC3K subset.
#'
#' @return A Seurat object from [load_demo_pbmc()].
#' @export
load_pbmc3k_demo <- function() {
  warning(
    "load_pbmc3k_demo() is a legacy alias. Returning load_demo_pbmc().",
    call. = FALSE
  )
  load_demo_pbmc()
}
