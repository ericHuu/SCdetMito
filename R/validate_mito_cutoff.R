# SCdetMito
# Author: Silu Hu
# Contact: husilu0902@gmail.com
# Version: 1.4.5
# Last updated: 2026-09-12

#' Validate an existing mitochondrial cutoff against basic QC metrics
#'
#' @description
#' `validate_mito_cutoff()` provides a deliberately secondary check of an
#' already selected mitochondrial cutoff. It does not estimate, move, or
#' replace that cutoff. Instead, it asks whether cells above the cutoff also
#' show lower detected-feature and count summaries, which is a common signature
#' of damaged or low-complexity cells. Discordance is retained as evidence for
#' review rather than used to silently alter the threshold.
#'
#' @param seurat_obj A Seurat object with mitochondrial and QC metadata.
#' @param cutoff A single mitochondrial cutoff, expressed as a fraction by
#'   default. Percentage-style values can be supplied with
#'   `cutoff_scale = "percent"`.
#' @param mito_col Metadata column containing mitochondrial fractions.
#' @param sample_col Optional metadata column for sample-stratified validation.
#'   The returned table always includes an overall row.
#' @param nFeature_col Metadata column containing detected feature counts.
#' @param nCount_col Metadata column containing total counts or UMIs.
#' @param min_cells_per_group Minimum number of cells required both at or below
#'   and above the cutoff for a conclusive comparison.
#' @param min_relative_drop Minimum relative decrease required in both
#'   `nFeature_col` and `nCount_col` to label the evidence `"supportive"`.
#' @param cutoff_scale Interpretation of `cutoff`: `"fraction"`, `"percent"`,
#'   or automatic scale detection.
#'
#' @return A data frame containing overall and optional sample-level cell
#'   counts, metric medians, relative changes, and an evidence classification.
#' @export
validate_mito_cutoff <- function(seurat_obj,
                                 cutoff,
                                 mito_col = "mitoRatio",
                                 sample_col = NULL,
                                 nFeature_col = "nFeature_RNA",
                                 nCount_col = "nCount_RNA",
                                 min_cells_per_group = 20,
                                 min_relative_drop = 0.10,
                                 cutoff_scale = c("fraction", "percent", "auto")) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.", call. = FALSE)
  }
  cutoff_scale <- match.arg(cutoff_scale)
  cutoff <- normalize_mito_cutoff_value(cutoff, "cutoff", scale = cutoff_scale)
  validate_numeric_scalar(min_cells_per_group, "min_cells_per_group", lower = 1)
  validate_numeric_scalar(min_relative_drop, "min_relative_drop", lower = 0, upper = 1)

  metadata <- seurat_obj@meta.data
  required <- c(mito_col, if (!is.null(sample_col)) sample_col)
  missing_required <- setdiff(required, colnames(metadata))
  if (length(missing_required)) {
    stop(
      "Missing required metadata column(s): ",
      paste(missing_required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  available_metrics <- intersect(c(nFeature_col, nCount_col), colnames(metadata))
  if (!length(available_metrics)) {
    stop(
      "Neither QC metric column is available: '", nFeature_col, "' and '", nCount_col, "'.",
      call. = FALSE
    )
  }

  mito <- as.numeric(metadata[[mito_col]])
  if (any(!is.finite(mito))) {
    stop("The mitochondrial ratio column must contain only finite values.", call. = FALSE)
  }

  groups <- list(overall = rep(TRUE, nrow(metadata)))
  if (!is.null(sample_col)) {
    sample_values <- as.character(metadata[[sample_col]])
    if (anyNA(sample_values)) {
      stop("The sample column cannot contain missing values.", call. = FALSE)
    }
    for (sample_id in unique(sample_values)) {
      groups[[paste0("sample:", sample_id)]] <- sample_values == sample_id
    }
  }

  relative_change <- function(above, below) {
    if (!is.finite(above) || !is.finite(below)) return(NA_real_)
    (above - below) / max(abs(below), .Machine$double.eps)
  }
  metric_median <- function(values, idx) {
    values <- suppressWarnings(as.numeric(values[idx]))
    values <- values[is.finite(values)]
    if (length(values)) stats::median(values) else NA_real_
  }

  rows <- lapply(names(groups), function(group_name) {
    in_group <- groups[[group_name]]
    below <- in_group & mito <= cutoff
    above <- in_group & mito > cutoff
    n_below <- sum(below)
    n_above <- sum(above)

    feature_below <- if (nFeature_col %in% available_metrics) {
      metric_median(metadata[[nFeature_col]], below)
    } else NA_real_
    feature_above <- if (nFeature_col %in% available_metrics) {
      metric_median(metadata[[nFeature_col]], above)
    } else NA_real_
    count_below <- if (nCount_col %in% available_metrics) {
      metric_median(metadata[[nCount_col]], below)
    } else NA_real_
    count_above <- if (nCount_col %in% available_metrics) {
      metric_median(metadata[[nCount_col]], above)
    } else NA_real_
    feature_change <- relative_change(feature_above, feature_below)
    count_change <- relative_change(count_above, count_below)

    complexity_below <- NA_real_
    complexity_above <- NA_real_
    if (all(c(nFeature_col, nCount_col) %in% available_metrics)) {
      feature_values <- suppressWarnings(as.numeric(metadata[[nFeature_col]]))
      count_values <- suppressWarnings(as.numeric(metadata[[nCount_col]]))
      complexity <- log10(pmax(feature_values, 0) + 1) /
        pmax(log10(pmax(count_values, 0) + 1), .Machine$double.eps)
      complexity_below <- metric_median(complexity, below)
      complexity_above <- metric_median(complexity, above)
    }

    enough_cells <- n_below >= min_cells_per_group && n_above >= min_cells_per_group
    has_both_primary_metrics <- all(is.finite(c(feature_change, count_change)))
    evidence_status <- if (!enough_cells || !has_both_primary_metrics) {
      "unavailable"
    } else if (feature_change <= -min_relative_drop && count_change <= -min_relative_drop) {
      "supportive"
    } else if (feature_change >= min_relative_drop && count_change >= min_relative_drop) {
      "discordant"
    } else {
      "inconclusive"
    }
    interpretation <- switch(
      evidence_status,
      supportive = "Cells above the supplied cutoff have lower median feature and count values.",
      discordant = "Cells above the supplied cutoff have higher median feature and count values.",
      inconclusive = "Feature and count changes are mixed or smaller than the predefined effect threshold.",
      unavailable = "Too few cells or too few available QC metrics for a conclusive comparison."
    )

    data.frame(
      scope = if (identical(group_name, "overall")) "overall" else "sample",
      sample = if (identical(group_name, "overall")) NA_character_ else sub("^sample:", "", group_name),
      cutoff = cutoff,
      n_cells = sum(in_group),
      n_at_or_below = n_below,
      n_above = n_above,
      fraction_at_or_below = if (sum(in_group)) n_below / sum(in_group) else NA_real_,
      nFeature_median_at_or_below = feature_below,
      nFeature_median_above = feature_above,
      nFeature_relative_change = feature_change,
      nCount_median_at_or_below = count_below,
      nCount_median_above = count_above,
      nCount_relative_change = count_change,
      complexity_median_at_or_below = complexity_below,
      complexity_median_above = complexity_above,
      evidence_status = evidence_status,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
