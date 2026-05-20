#' load_demo_seurat: load the bundled multi-group Seurat demo object
#'
#' @return A Seurat object stored in `inst/extdata/demo_multigroup_seurat.rds`.
#'   The demo contains six samples across three biological groups and includes a
#'   precomputed `mitoRatio` metadata column.
#' @export
load_demo_seurat <- function() {
  demo_path <- system.file("extdata", "demo_multigroup_seurat.rds", package = "SCdetMito")
  if (!nzchar(demo_path)) {
    stop("The bundled demo Seurat object could not be located.", call. = FALSE)
  }

  readRDS(demo_path)
}

#' load_pbmc3k_demo: load the bundled lightweight public PBMC3K object
#'
#' @description
#' `load_pbmc3k_demo()` loads a small subset derived from the public 10x
#' Genomics PBMC3K dataset. It is intended for examples and package smoke tests
#' without making the installed package large.
#'
#' @return A Seurat object stored in `inst/extdata/pbmc3k_lightweight_seurat.rds`.
#'   The object contains sampled PBMC3K cells, a small expression feature set,
#'   pseudo-sample labels for multi-sample QC examples, and a precomputed
#'   `mitoRatio` metadata column.
#' @export
load_pbmc3k_demo <- function() {
  demo_path <- system.file("extdata", "pbmc3k_lightweight_seurat.rds", package = "SCdetMito")
  if (!nzchar(demo_path)) {
    stop(
      "The bundled lightweight PBMC3K Seurat object could not be located.",
      call. = FALSE
    )
  }

  readRDS(demo_path)
}
