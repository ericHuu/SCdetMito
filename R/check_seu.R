#' check_seu: validate that a metadata column exists in a Seurat object
#'
#' @param seurat_obj A Seurat object.
#' @param check Metadata column name to validate. Defaults to `"mitoRatio"`.
#' @param normalize_fraction Whether to rescale percentage-like values to
#'   fractions in the 0-1 range. Defaults to `FALSE`.
#' @param must_be_numeric Whether the target metadata column must be numeric.
#'   Defaults to `FALSE`.
#'
#' @return The input Seurat object, potentially with a normalized metadata
#'   column.
#' @export
#'
#' @examples
#' # DO NOT RUN
#' # library(Seurat)
#' # seu <- CreateSeuratObject(matrix(rpois(2000, 5), nrow = 100))
#' # seu$mitoRatio <- runif(ncol(seu))
#' # check_seu(seu, "mitoRatio", normalize_fraction = TRUE)
#' # DO NOT RUN
check_seu <- function(seurat_obj,
                      check = "mitoRatio",
                      normalize_fraction = FALSE,
                      must_be_numeric = FALSE) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.", call. = FALSE)
  }

  if (!check %in% colnames(seurat_obj@meta.data)) {
    stop(
      "The specified column '", check, "' is not present in the Seurat object metadata.",
      call. = FALSE
    )
  }

  values <- seurat_obj@meta.data[[check]]

  if (must_be_numeric && !is.numeric(values)) {
    stop(
      "The specified column '", check, "' must be numeric.",
      call. = FALSE
    )
  }

  if (normalize_fraction && is.numeric(values)) {
    finite_values <- values[is.finite(values)]
    if (length(finite_values) && max(finite_values, na.rm = TRUE) > 1) {
      seurat_obj@meta.data[[check]] <- values / 100
    }
  }

  seurat_obj
}
