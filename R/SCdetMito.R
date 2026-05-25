# SCdetMito
# Author: Silu Hu
# Contact: husilu0902@gmail.com
# Version: 1.4.3
# Last updated: 2026-05-23

#' SCdetMito: reference-aware retention-loss mitochondrial QC
#'
#' @description
#' `SCdetMito()` estimates sample-aware mitochondrial ratio cutoffs from
#' retained-cell profiles and reports multiple interpretable QC decision
#' candidates. The function is kept under its original public name for backward
#' compatibility, but the statistical framing is reference-aware
#' retention-loss enrichment-based adaptive mitochondrial QC rather than a
#' breakpoint model.
#'
#' @details
#' For each sample `s` and candidate mitochondrial cutoff `c`, the retained-cell
#' profile is calculated as `R_s(c) = sum(m_i < c)`, where `m_i` is the
#' mitochondrial ratio of cell `i`. Adjacent candidate cutoffs define an
#' interval-specific cell loss, `L_s(c_j) = R_s(c_{j-1}) - R_s(c_j)`. Candidate
#' intervals are first filtered by absolute cell loss, relative cell loss,
#' retained cell count, and retained fraction. Retention-loss enrichment is then
#' evaluated within each sample, p-values are adjusted per sample, and
#' significant intervals are summarized into sample-level cutoff candidates:
#' the first significant high-to-low boundary, the largest-drop boundary, the
#' selected cutoff under the requested decision mode, and a package-level
#' `recommended_cutoff` derived by the default reference-aware recommendation
#' policy. When species/tissue or a user reference is supplied,
#' literature-informed reference cutoffs are added for interpretation and
#' warning generation. A sample-supported global cutoff is returned for
#' multi-sample inputs.
#' Confidence labels describe QC decision confidence only: high/moderate labels
#' require a statistically detected, non-fallback interval, whereas low labels
#' mark fallback-derived or not-detected cutoffs.
#'
#' @param seurat_obj A Seurat object.
#' @param mito_col,mitoRatio Metadata column storing mitochondrial ratios.
#'   `mito_col` is the preferred name; `mitoRatio` is retained for backward
#'   compatibility.
#' @param sample_col,by Metadata column identifying samples. `sample_col` is the
#'   preferred name; `by` is retained for backward compatibility.
#' @param bin_width Candidate cutoff grid step size.
#' @param min_cut,max_cut Minimum and maximum mitochondrial ratio cutoffs
#'   considered.
#' @param write_tables,table_out Whether to export retained-cell and interval
#'   loss tables. `write_tables` is preferred.
#' @param write_plots,plot Whether to export manuscript-ready PDF plots.
#'   `write_plots` is preferred.
#' @param min_drop_cells,diff_num Minimum absolute interval-specific cell loss.
#'   `min_drop_cells` is preferred; `diff_num` is retained for compatibility.
#' @param min_drop_fraction Minimum interval loss fraction relative to the
#'   sample cell count.
#' @param min_cells_after Minimum retained cells after the candidate cutoff.
#' @param min_retention_after Minimum retained fraction after the candidate
#'   cutoff.
#' @param min_sample_support,min_group_support Minimum fraction of samples that
#'   must support a global cutoff. `min_sample_support` is preferred.
#' @param loss_test Retention-loss enrichment test. Defaults to
#'   `"mad_zscore"` as the primary robust FDR-compatible detector. Other
#'   options are `"empirical_tail"`, `"poisson_tail"`, `"zscore"`, and
#'   `"threshold_only"`.
#' @param p_adj_method,p_adjust_method Multiple-testing correction method passed
#'   to [stats::p.adjust()]. Use `"none"` to leave p-values unadjusted.
#'   `p_adj_method` is preferred.
#' @param alpha Adjusted p-value threshold. Defaults to `0.05`.
#' @param sample_cutoff_method Rule for selecting one sample-level cutoff from
#'   significant intervals. Defaults to `"largest_drop"` for backward
#'   compatibility. `"reference_guided"` is recommended when a literature
#'   reference is supplied and the user wants the first significant high-to-low
#'   boundary prioritized over the largest-drop boundary.
#' @param species,tissue Optional species and tissue names used to look up
#'   literature-informed reference cutoffs. These references are decision
#'   priors, not hard filtering rules.
#' @param reference_cutoff Optional user-provided mitochondrial fraction
#'   reference cutoff. Overrides `species`/`tissue` lookup when supplied.
#' @param reference_table Optional custom reference table with columns
#'   `species`, `tissue`, and `reference_cutoff`.
#' @param reference_warning Whether to add warning flags/messages when selected
#'   data-driven cutoffs deviate strongly from available reference values or
#'   retain very few cells.
#' @param auto_add_mito Whether to calculate `mito_col` with
#'   [ensure_mito_ratio()] when the column is absent. Defaults to `TRUE`.
#' @param mito_features Optional complete mitochondrial feature vector forwarded
#'   to [ensure_mito_ratio()].
#' @param mito_pattern Optional mitochondrial feature regex forwarded to
#'   [ensure_mito_ratio()].
#' @param mito_assay Optional assay used for mitochondrial ratio calculation.
#' @param recompute_mito Whether to recompute `mito_col` even when it already
#'   exists. Defaults to `FALSE`, which preserves an existing validated
#'   mitochondrial fraction column.
#' @param output_dir Directory for exported outputs.
#' @param fallback_method Fallback used only when no significant interval is
#'   found for a sample. Options are `"quantile"`, `"max_cut"`, and `"none"`.
#' @param fallback_quantile Quantile used when `fallback_method = "quantile"`.
#' @param return_details Whether to return detailed tables. Defaults to `FALSE`.
#' @param ... Additional arguments ignored for backward compatibility.
#'
#' @return A numeric sample-supported global cutoff, or a list with detailed
#'   retained-cell profiles, interval-loss tables, significant intervals, and
#'   settings when `return_details = TRUE`. In `sample_cutoff_summary`,
#'   `selected_cutoff` is determined by the user-selected
#'   `sample_cutoff_method`, `detected_cutoff` is retained as a
#'   backward-compatible alias of `selected_cutoff`, and `recommended_cutoff`
#'   is SCdetMito's reference-aware recommendation reported regardless of the
#'   selected mode. The
#'   `warning_level` column is human-readable; logical warning columns such as
#'   `has_reference_deviation`, `has_low_retention_warning`, and
#'   `has_fallback_warning` are intended for programmatic use. The
#'   `change_points` element is retained as a legacy alias of
#'   `interval_loss_table`.
#' @importFrom rlang .data
#' @export
SCdetMito <- function(seurat_obj,
                      mito_col = NULL,
                      sample_col = NULL,
                      mitoRatio = "mitoRatio",
                      by = "sample",
                      bin_width = 0.01,
                      min_cut = 0.01,
                      max_cut = 1,
                      write_tables = NULL,
                      write_plots = NULL,
                      table_out = FALSE,
                      plot = TRUE,
                      min_drop_cells = NULL,
                      diff_num = 5,
                      min_drop_fraction = 0.01,
                      min_cells_after = 20,
                      min_retention_after = 0.1,
                      min_sample_support = NULL,
                      min_group_support = 0.5,
                      loss_test = c("mad_zscore", "empirical_tail", "poisson_tail", "zscore", "threshold_only"),
                      p_adj_method = NULL,
                      p_adjust_method = "BH",
                      alpha = 0.05,
                      sample_cutoff_method = c(
                        "largest_drop",
                        "first_significant_high",
                        "first_significant_low",
                        "reference_guided",
                        "largest_drop_with_reference_guard",
                        "max_significant",
                        "median_significant",
                        "min_significant"
                      ),
                      species = NULL,
                      tissue = NULL,
                      reference_cutoff = NULL,
                      reference_table = NULL,
                      reference_warning = TRUE,
                      auto_add_mito = TRUE,
                      mito_features = NULL,
                      mito_pattern = NULL,
                      mito_assay = NULL,
                      recompute_mito = FALSE,
                      output_dir = ".",
                      fallback_method = c("quantile", "max_cut", "none"),
                      fallback_quantile = 0.9,
                      return_details = FALSE,
                      ...) {
  message("Estimating retention-loss enrichment-based mitochondrial cutoffs...")

  mito_col <- mito_col %||% mitoRatio
  sample_col <- sample_col %||% by
  write_tables <- write_tables %||% table_out
  write_plots <- write_plots %||% plot
  min_drop_cells <- min_drop_cells %||% diff_num
  min_sample_support <- min_sample_support %||% min_group_support
  p_adj_method <- p_adj_method %||% p_adjust_method

  seurat_obj <- check_seu(seurat_obj, sample_col)
  if (!mito_col %in% colnames(seurat_obj@meta.data) && !isTRUE(auto_add_mito)) {
    stop(
      "Mitochondrial column '",
      mito_col,
      "' is missing and auto_add_mito = FALSE. Provide a valid mito_col or set auto_add_mito = TRUE.",
      call. = FALSE
    )
  }
  seurat_obj <- ensure_mito_ratio(
    object = seurat_obj,
    mito_col = mito_col,
    species = species,
    assay = mito_assay,
    mito_features = mito_features,
    mito_pattern = mito_pattern,
    overwrite = isTRUE(recompute_mito),
    convert_percent = TRUE,
    verbose = FALSE
  )
  output_dir <- ensure_output_dir(output_dir)
  fallback_method <- match.arg(fallback_method)
  loss_test <- match.arg(loss_test)
  sample_cutoff_method <- match.arg(sample_cutoff_method)
  reference_warning <- isTRUE(reference_warning)
  if (identical(loss_test, "empirical_tail") && !identical(p_adj_method, "none")) {
    warning(
      "empirical_tail p-values are rank-based and highly discrete; multiple-testing adjustment may be overly conservative and can increase fallback frequency.",
      call. = FALSE
    )
  }

  validate_numeric_scalar(bin_width, "bin_width", lower = .Machine$double.eps)
  if (!is.numeric(min_cut) || !is.numeric(max_cut) || min_cut >= max_cut) {
    stop("'min_cut' must be smaller than 'max_cut'.", call. = FALSE)
  }
  validate_numeric_scalar(min_drop_cells, "min_drop_cells", lower = 0)
  validate_numeric_scalar(min_drop_fraction, "min_drop_fraction", lower = 0, upper = 1)
  validate_numeric_scalar(min_cells_after, "min_cells_after", lower = 1)
  validate_numeric_scalar(min_retention_after, "min_retention_after", lower = 0, upper = 1)
  validate_numeric_scalar(min_sample_support, "min_sample_support", lower = .Machine$double.eps, upper = 1)
  validate_numeric_scalar(alpha, "alpha", lower = .Machine$double.eps, upper = 1)
  validate_numeric_scalar(fallback_quantile, "fallback_quantile", lower = .Machine$double.eps, upper = 1)
  reference_info <- .get_reference_cutoff(
    species = species,
    tissue = tissue,
    default_reference = reference_cutoff,
    reference_table = reference_table,
    verbose = TRUE
  )
  if (!is.character(p_adj_method) ||
    length(p_adj_method) != 1 ||
    is.na(p_adj_method) ||
    !p_adj_method %in% c("none", stats::p.adjust.methods)) {
    stop(
      "'p_adj_method' must be one of: ",
      paste(c("none", stats::p.adjust.methods), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  metadata <- seurat_obj@meta.data
  samples <- as.character(metadata[[sample_col]])
  mito_values <- as.numeric(metadata[[mito_col]])
  if (anyNA(samples)) {
    stop("The sample column cannot contain missing values.", call. = FALSE)
  }
  if (anyNA(mito_values)) {
    stop("The mitochondrial ratio column cannot contain missing values.", call. = FALSE)
  }

  sample_levels <- unique(samples)
  n_cells_by_sample <- as.integer(table(factor(samples, levels = sample_levels)))
  names(n_cells_by_sample) <- sample_levels

  snap_cutoff <- function(value) {
    if (!is.finite(value)) {
      return(NA_real_)
    }
    snapped <- round(value / bin_width) * bin_width
    snapped <- max(min_cut, min(max_cut, snapped))
    round(snapped, 6)
  }

  cutoff_grid <- round(seq(from = max_cut, to = min_cut, by = -bin_width), 6)
  cutoff_grid <- cutoff_grid[cutoff_grid >= min_cut & cutoff_grid <= max_cut]
  if (!length(cutoff_grid) || length(cutoff_grid) < 2) {
    stop("The cutoff grid must contain at least two candidate cutoffs.", call. = FALSE)
  }

  retained_cell_profile <- do.call(
    rbind,
    lapply(sample_levels, function(sample_id) {
      values <- mito_values[samples == sample_id]
      retained <- vapply(cutoff_grid, function(cutoff) {
        sum(values < cutoff, na.rm = TRUE)
      }, integer(1))
      data.frame(
        sample = sample_id,
        cutoff = cutoff_grid,
        n_cells = length(values),
        retained_cells = retained,
        retention_fraction = retained / length(values),
        stringsAsFactors = FALSE
      )
    })
  )

  retained_cell_profile_legacy <- reshape2::dcast(
    retained_cell_profile,
    sample + n_cells ~ cutoff,
    value.var = "retained_cells"
  )
  names(retained_cell_profile_legacy)[names(retained_cell_profile_legacy) == "sample"] <- sample_col
  names(retained_cell_profile_legacy)[names(retained_cell_profile_legacy) == "n_cells"] <- "raw_counts"

  interval_loss_table <- do.call(
    rbind,
    lapply(sample_levels, function(sample_id) {
      profile <- retained_cell_profile[retained_cell_profile$sample == sample_id, , drop = FALSE]
      profile <- profile[order(-profile$cutoff), , drop = FALSE]
      if (nrow(profile) < 2) {
        return(data.frame())
      }

      interval_loss <- pmax(profile$retained_cells[-nrow(profile)] - profile$retained_cells[-1], 0)
      candidate_cutoffs <- profile$cutoff[-1]
      retained_after <- profile$retained_cells[-1]
      retention_after <- profile$retention_fraction[-1]
      n_cells <- profile$n_cells[[1]]

      rows <- lapply(seq_along(interval_loss), function(i) {
        background_losses <- interval_loss[-i]
        test_result <- run_loss_enrichment_test(
          loss_value = interval_loss[[i]],
          background_losses = background_losses,
          method = loss_test
        )
        interval_loss_fraction <- interval_loss[[i]] / n_cells
        passes_min_drop_cells <- interval_loss[[i]] >= min_drop_cells
        passes_min_drop_fraction <- interval_loss_fraction >= min_drop_fraction
        passes_min_cells_after <- retained_after[[i]] >= min_cells_after
        passes_min_retention_after <- retention_after[[i]] >= min_retention_after
        data.frame(
          sample = sample_id,
          interval_index = i,
          previous_cutoff = profile$cutoff[[i]],
          cutoff = candidate_cutoffs[[i]],
          n_cells = n_cells,
          retained_cells = retained_after[[i]],
          retention_fraction = retention_after[[i]],
          interval_loss = interval_loss[[i]],
          interval_loss_fraction = interval_loss_fraction,
          passes_min_drop_cells = passes_min_drop_cells,
          passes_min_drop_fraction = passes_min_drop_fraction,
          passes_min_cells_after = passes_min_cells_after,
          passes_min_retention_after = passes_min_retention_after,
          passes_candidate_filters = passes_min_drop_cells &&
            passes_min_drop_fraction &&
            passes_min_cells_after &&
            passes_min_retention_after,
          p_value = test_result$p_value,
          p_adj = NA_real_,
          loss_score = test_result$loss_score,
          test_method = test_result$test_method,
          baseline_loss = test_result$baseline_loss,
          is_significant = FALSE,
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, rows)
    })
  )

  if (!nrow(interval_loss_table)) {
    interval_loss_table <- data.frame(
      sample = character(),
      interval_index = integer(),
      previous_cutoff = numeric(),
      cutoff = numeric(),
      n_cells = integer(),
      retained_cells = integer(),
      retention_fraction = numeric(),
      interval_loss = numeric(),
      interval_loss_fraction = numeric(),
      passes_min_drop_cells = logical(),
      passes_min_drop_fraction = logical(),
      passes_min_cells_after = logical(),
      passes_min_retention_after = logical(),
      passes_candidate_filters = logical(),
      p_value = numeric(),
      p_adj = numeric(),
      loss_score = numeric(),
      test_method = character(),
      baseline_loss = numeric(),
      is_significant = logical(),
      stringsAsFactors = FALSE
    )
  }

  interval_loss_table$p_adj_method <- p_adj_method
  if (nrow(interval_loss_table)) {
    for (sample_id in sample_levels) {
      idx <- which(interval_loss_table$sample == sample_id)
      if (identical(loss_test, "threshold_only")) {
        interval_loss_table$is_significant[idx] <- interval_loss_table$passes_candidate_filters[idx]
      } else if (identical(p_adj_method, "none")) {
        interval_loss_table$p_adj[idx] <- interval_loss_table$p_value[idx]
        interval_loss_table$is_significant[idx] <-
          interval_loss_table$passes_candidate_filters[idx] &
            interval_loss_table$p_adj[idx] <= alpha
      } else {
        interval_loss_table$p_adj[idx] <- stats::p.adjust(
          interval_loss_table$p_value[idx],
          method = p_adj_method
        )
        interval_loss_table$is_significant[idx] <-
          interval_loss_table$passes_candidate_filters[idx] &
            interval_loss_table$p_adj[idx] <= alpha
      }
    }
  }
  interval_loss_table$p_adjusted <- interval_loss_table$p_adj
  interval_loss_table$p_fdr <- interval_loss_table$p_adj
  interval_loss_table$p_adjust_method <- interval_loss_table$p_adj_method
  interval_loss_table$drop_fraction <- interval_loss_table$interval_loss_fraction
  interval_loss_table$delta_prev <- interval_loss_table$interval_loss

  significant_intervals <- interval_loss_table[
    interval_loss_table$is_significant %in% TRUE,
    ,
    drop = FALSE
  ]

  fallback_cutoff_for_sample <- function(sample_id) {
    if (identical(fallback_method, "none")) {
      return(NA_real_)
    }
    if (identical(fallback_method, "max_cut")) {
      return(max_cut)
    }
    values <- mito_values[samples == sample_id]
    snap_cutoff(as.numeric(stats::quantile(
      values,
      probs = fallback_quantile,
      na.rm = TRUE,
      names = FALSE,
      type = 8
    )))
  }

  source_for_fallback <- function() {
    if (identical(fallback_method, "quantile")) {
      paste0("quantile_fallback_", format(fallback_quantile, nsmall = 2, trim = TRUE))
    } else if (identical(fallback_method, "max_cut")) {
      "max_cut_fallback"
    } else {
      "not_detected"
    }
  }

  sample_cutoff_summary <- do.call(
    rbind,
    lapply(sample_levels, function(sample_id) {
      sample_table <- interval_loss_table[interval_loss_table$sample == sample_id, , drop = FALSE]
      sample_sig <- sample_table[sample_table$is_significant %in% TRUE, , drop = FALSE]
      fallback_cutoff <- fallback_cutoff_for_sample(sample_id)
      fallback_used <- !nrow(sample_sig) && !identical(fallback_method, "none")

      first_significant_cutoff_high <- if (nrow(sample_sig)) {
        snap_cutoff(max(sample_sig$cutoff, na.rm = TRUE))
      } else {
        NA_real_
      }
      first_significant_cutoff_low <- if (nrow(sample_sig)) {
        snap_cutoff(min(sample_sig$cutoff, na.rm = TRUE))
      } else {
        NA_real_
      }
      largest_drop_cutoff <- if (nrow(sample_sig)) {
        ranked_sig <- sample_sig[order(-sample_sig$interval_loss, -sample_sig$cutoff), , drop = FALSE]
        snap_cutoff(ranked_sig$cutoff[[1]])
      } else {
        NA_real_
      }

      if (nrow(sample_sig)) {
        selected_cutoff <- switch(sample_cutoff_method,
          largest_drop = largest_drop_cutoff,
          largest_drop_with_reference_guard = largest_drop_cutoff,
          first_significant_high = first_significant_cutoff_high,
          reference_guided = if (is.finite(first_significant_cutoff_high)) {
            first_significant_cutoff_high
          } else {
            largest_drop_cutoff
          },
          first_significant_low = first_significant_cutoff_low,
          max_significant = snap_cutoff(max(sample_sig$cutoff, na.rm = TRUE)),
          median_significant = snap_cutoff(stats::median(unique(sample_sig$cutoff), na.rm = TRUE)),
          min_significant = snap_cutoff(min(sample_sig$cutoff, na.rm = TRUE))
        )
        cutoff_source <- switch(sample_cutoff_method,
          largest_drop = "significant_largest_drop",
          largest_drop_with_reference_guard = "significant_largest_drop_with_reference_guard",
          first_significant_high = "first_significant_high",
          reference_guided = if (identical(selected_cutoff, first_significant_cutoff_high)) {
            "reference_guided_first_significant_high"
          } else {
            "reference_guided_largest_drop"
          },
          first_significant_low = "first_significant_low",
          max_significant = "significant_max_significant",
          median_significant = "significant_median_significant",
          min_significant = "significant_min_significant"
        )
        selected_cutoff_interpretation <- switch(sample_cutoff_method,
          largest_drop = "Largest interval-specific cell loss among significant intervals.",
          largest_drop_with_reference_guard = "Largest interval-specific cell loss among significant intervals, with reference and retention warning fields reported.",
          first_significant_high = "First significant high-to-low retention-loss boundary.",
          reference_guided = "Reference-aware decision mode prioritizing the first significant high-to-low boundary when available.",
          first_significant_low = "First significant low-to-high retention-loss boundary.",
          max_significant = "Largest significant cutoff.",
          median_significant = "Median significant cutoff snapped to the cutoff grid.",
          min_significant = "Smallest significant cutoff."
        )
      } else {
        selected_cutoff <- fallback_cutoff
        cutoff_source <- source_for_fallback()
        selected_cutoff_interpretation <- if (is.finite(selected_cutoff)) {
          "Fallback-derived cutoff because no significant interval-specific cell-loss boundary was detected."
        } else {
          "No cutoff detected and fallback_method = 'none'."
        }
      }
      detected_cutoff <- selected_cutoff

      find_profile_at_cutoff <- function(cutoff_value) {
        profile_at_cutoff <- retained_cell_profile[
          retained_cell_profile$sample == sample_id &
            retained_cell_profile$cutoff == cutoff_value,
          ,
          drop = FALSE
        ]
        if (!nrow(profile_at_cutoff) && is.finite(cutoff_value)) {
          profile_at_cutoff <- retained_cell_profile[
            retained_cell_profile$sample == sample_id &
              retained_cell_profile$cutoff == snap_cutoff(cutoff_value),
            ,
            drop = FALSE
          ]
        }
        profile_at_cutoff
      }
      profile_at_cutoff <- find_profile_at_cutoff(detected_cutoff)
      recommendation_preview <- build_recommendation_fields(
        first_significant_cutoff_high = first_significant_cutoff_high,
        largest_drop_cutoff = largest_drop_cutoff,
        fallback_cutoff = fallback_cutoff,
        fallback_used = fallback_used,
        fallback_method = fallback_method,
        fallback_quantile = fallback_quantile,
        reference_cutoff = reference_info$reference_cutoff,
        retention_fraction_at_recommended = NA_real_,
        reference_warning = reference_warning
      )
      profile_at_recommended <- find_profile_at_cutoff(recommendation_preview$recommended_cutoff)

      largest_interval_loss <- if (nrow(sample_table)) max(sample_table$interval_loss, na.rm = TRUE) else NA_real_
      largest_interval_loss_fraction <- if (nrow(sample_table)) {
        max(sample_table$interval_loss_fraction, na.rm = TRUE)
      } else {
        NA_real_
      }
      cutoff_confidence <- classify_cutoff_confidence(
        fallback_used = fallback_used,
        cutoff_source = cutoff_source,
        significant_interval_count = nrow(sample_sig),
        largest_interval_loss_fraction = largest_interval_loss_fraction
      )
      retained_cells_at_cutoff <- if (nrow(profile_at_cutoff)) profile_at_cutoff$retained_cells[[1]] else NA_real_
      retention_fraction_at_cutoff <- if (nrow(profile_at_cutoff)) profile_at_cutoff$retention_fraction[[1]] else NA_real_
      reference_ratio <- if (is.finite(selected_cutoff) &&
        is.finite(reference_info$reference_cutoff) &&
        reference_info$reference_cutoff > 0) {
          selected_cutoff / reference_info$reference_cutoff
      } else {
        NA_real_
      }
      reference_deviation_flag <- is.finite(reference_ratio) && reference_ratio >= 2
      recommended_retained_cells <- if (nrow(profile_at_recommended)) profile_at_recommended$retained_cells[[1]] else NA_real_
      recommended_retention_fraction <- if (nrow(profile_at_recommended)) profile_at_recommended$retention_fraction[[1]] else NA_real_
      recommendation_info <- build_recommendation_fields(
        first_significant_cutoff_high = first_significant_cutoff_high,
        largest_drop_cutoff = largest_drop_cutoff,
        fallback_cutoff = fallback_cutoff,
        fallback_used = fallback_used,
        fallback_method = fallback_method,
        fallback_quantile = fallback_quantile,
        reference_cutoff = reference_info$reference_cutoff,
        retention_fraction_at_recommended = recommended_retention_fraction,
        reference_warning = reference_warning
      )
      warning_info <- build_reference_warning(
        selected_cutoff = selected_cutoff,
        reference_cutoff = reference_info$reference_cutoff,
        retention_fraction_at_cutoff = retention_fraction_at_cutoff,
        fallback_used = fallback_used,
        sample_cutoff_method = sample_cutoff_method,
        first_significant_cutoff_high = first_significant_cutoff_high,
        largest_drop_cutoff = largest_drop_cutoff,
        reference_warning = reference_warning
      )

      data.frame(
        sample = sample_id,
        n_cells = n_cells_by_sample[[sample_id]],
        detected_cutoff = detected_cutoff,
        selected_cutoff = selected_cutoff,
        selected_cutoff_method = sample_cutoff_method,
        selected_cutoff_source = cutoff_source,
        selected_cutoff_interpretation = selected_cutoff_interpretation,
        recommended_cutoff = recommendation_info$recommended_cutoff,
        recommended_method = recommendation_info$recommended_method,
        recommended_reason = recommendation_info$recommended_reason,
        recommendation_level = recommendation_info$recommendation_level,
        recommendation_warning = recommendation_info$recommendation_warning,
        recommendation_source = recommendation_info$recommendation_source,
        cutoff_source = cutoff_source,
        fallback_used = fallback_used,
        fallback_method = if (fallback_used) fallback_method else NA_character_,
        fallback_quantile = if (fallback_used && identical(fallback_method, "quantile")) fallback_quantile else NA_real_,
        significant_interval_count = nrow(sample_sig),
        significant_cutoff_count = nrow(sample_sig),
        cutoff_confidence = cutoff_confidence,
        reference_cutoff = reference_info$reference_cutoff,
        reference_species = reference_info$reference_species,
        reference_tissue = reference_info$reference_tissue,
        reference_source = reference_info$reference_source,
        first_significant_cutoff_high = first_significant_cutoff_high,
        first_significant_cutoff_low = first_significant_cutoff_low,
        largest_drop_cutoff = largest_drop_cutoff,
        reference_ratio = reference_ratio,
        reference_deviation_flag = reference_deviation_flag,
        warning_level = warning_info$warning_level,
        warning_message = warning_info$warning_message,
        has_reference_deviation = warning_info$has_reference_deviation,
        has_low_retention_warning = warning_info$has_low_retention_warning,
        has_fallback_warning = warning_info$has_fallback_warning,
        has_high_reference_deviation = warning_info$has_high_reference_deviation,
        has_moderate_reference_deviation = warning_info$has_moderate_reference_deviation,
        largest_interval_loss = largest_interval_loss,
        largest_interval_loss_fraction = largest_interval_loss_fraction,
        retained_cells_at_cutoff = retained_cells_at_cutoff,
        retention_fraction_at_cutoff = retention_fraction_at_cutoff,
        recommended_retained_cells = recommended_retained_cells,
        recommended_retention_fraction = recommended_retention_fraction,
        loss_test = loss_test,
        sample_cutoff_method = sample_cutoff_method,
        stringsAsFactors = FALSE
      )
    })
  )

  detected_cutoffs <- sample_cutoff_summary$detected_cutoff
  finite_cutoffs <- detected_cutoffs[is.finite(detected_cutoffs)]
  if (length(finite_cutoffs)) {
    support <- resolve_supported_cutoff(
      detected_cutoffs,
      min_support_fraction = min_sample_support,
      n_total_samples = length(detected_cutoffs)
    )
    global_cutoff <- support$global_cutoff
  } else {
    support <- data.frame(
      global_cutoff = NA_real_,
      support_fraction = 0,
      n_supporting_samples = 0L,
      n_total_samples = length(detected_cutoffs),
      required_support = max(1L, ceiling(length(detected_cutoffs) * min_sample_support)),
      strategy = "sample_supported_global_cutoff",
      stringsAsFactors = FALSE
    )
    global_cutoff <- NA_real_
  }
  sample_cutoff_summary$supports_final_cutoff <- is.finite(sample_cutoff_summary$detected_cutoff) &
    is.finite(global_cutoff) &
    sample_cutoff_summary$detected_cutoff >= global_cutoff
  recommended_cutoffs <- sample_cutoff_summary$recommended_cutoff
  finite_recommended_cutoffs <- recommended_cutoffs[is.finite(recommended_cutoffs)]
  recommended_support <- if (length(finite_recommended_cutoffs)) {
    resolve_supported_cutoff(
      recommended_cutoffs,
      min_support_fraction = min_sample_support,
      n_total_samples = length(recommended_cutoffs),
      strategy = "recommended_sample_supported_global_cutoff"
    )
  } else {
    data.frame(
      global_cutoff = NA_real_,
      cutoff = NA_real_,
      support_fraction = 0,
      n_supporting_samples = 0L,
      support_count = 0L,
      n_total_samples = length(recommended_cutoffs),
      sample_count = length(recommended_cutoffs),
      required_support = max(1L, ceiling(length(recommended_cutoffs) * min_sample_support)),
      strategy = "recommended_sample_supported_global_cutoff",
      group = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  sample_cutoff_summary$supports_recommended_final_cutoff <- is.finite(sample_cutoff_summary$recommended_cutoff) &
    is.finite(recommended_support$global_cutoff) &
    sample_cutoff_summary$recommended_cutoff >= recommended_support$global_cutoff

  if (is.finite(global_cutoff)) {
    message("Sample-supported global mito cutoff: ", global_cutoff)
  } else {
    warning("No sample-supported global mitochondrial cutoff was detected.", call. = FALSE)
  }

  if (write_tables) {
    utils::write.csv(
      retained_cell_profile,
      file = file.path(output_dir, "SCdetMito_retained_cell_profile.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      interval_loss_table,
      file = file.path(output_dir, "SCdetMito_interval_loss_table.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      significant_intervals,
      file = file.path(output_dir, "SCdetMito_significant_intervals.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      sample_cutoff_summary,
      file = file.path(output_dir, "SCdetMito_sample_cutoff_summary.csv"),
      row.names = FALSE
    )
    reference_columns <- intersect(
      c(
        "sample", "reference_cutoff", "reference_species",
        "reference_tissue", "reference_source", "reference_ratio",
        "reference_deviation_flag", "warning_level", "warning_message",
        "has_reference_deviation", "has_high_reference_deviation",
        "has_moderate_reference_deviation"
      ),
      colnames(sample_cutoff_summary)
    )
    if (length(reference_columns) > 1L) {
      utils::write.csv(
        sample_cutoff_summary[, reference_columns, drop = FALSE],
        file = file.path(output_dir, "SCdetMito_reference_cutoff_summary.csv"),
        row.names = FALSE
      )
    }
    recommendation_columns <- intersect(
      c(
        "sample", "selected_cutoff", "detected_cutoff", "recommended_cutoff",
        "recommended_method", "recommendation_source", "recommendation_level",
        "recommended_reason", "recommendation_warning", "reference_cutoff",
        "first_significant_cutoff_high", "first_significant_cutoff_low",
        "largest_drop_cutoff", "retained_cells_at_cutoff",
        "retention_fraction_at_cutoff", "recommended_retained_cells",
        "recommended_retention_fraction", "warning_level", "warning_message"
      ),
      colnames(sample_cutoff_summary)
    )
    utils::write.csv(
      sample_cutoff_summary[, recommendation_columns, drop = FALSE],
      file = file.path(output_dir, "SCdetMito_recommendation_summary.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      retained_cell_profile_legacy,
      file = file.path(output_dir, "processed_temp_cell_num.csv"),
      quote = FALSE,
      row.names = FALSE
    )
    # Legacy interval-loss filename retained for backward compatibility.
    utils::write.csv(
      interval_loss_table,
      file = file.path(output_dir, "mito_change_point_results.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      sample_cutoff_summary,
      file = file.path(output_dir, "mito_sample_cutoff_summary.csv"),
      row.names = FALSE
    )
  }

  if (write_plots) {
    plot_scdetmito_outputs(
      retained_cell_profile = retained_cell_profile,
      interval_loss_table = interval_loss_table,
      sample_cutoff_summary = sample_cutoff_summary,
      output_dir = output_dir
    )
  }

  settings <- list(
    method_label = "retention-loss enrichment-based adaptive mitochondrial cutoff detection",
    retained_cell_profile_formula = "R_s(c) = sum(m_i < c)",
    interval_loss_formula = "L_s(c_j) = R_s(c_{j-1}) - R_s(c_j)",
    loss_test = loss_test,
    p_adj_method = p_adj_method,
    p_adjust_method = p_adj_method,
    alpha = alpha,
    sample_cutoff_method = sample_cutoff_method,
    recommendation_policy = "first_significant_high_then_largest_drop_then_fallback",
    detector_mode = classify_detector_mode(loss_test, p_adj_method),
    reference_cutoff = reference_info$reference_cutoff,
    reference_species = reference_info$reference_species,
        reference_tissue = reference_info$reference_tissue,
        reference_source = reference_info$reference_source,
        reference_note = reference_info$reference_note,
        reference_found = reference_info$reference_found,
        reference_warning = reference_warning,
        auto_add_mito = auto_add_mito,
        mito_features_provided = !is.null(mito_features),
        mito_pattern = if (is.null(mito_pattern)) NA_character_ else paste(mito_pattern, collapse = ";"),
        mito_assay = mito_assay %||% Seurat::DefaultAssay(seurat_obj),
        recompute_mito = isTRUE(recompute_mito),
        final_cutoff_rule = "sample_supported_global_cutoff",
    min_sample_support = min_sample_support,
    min_group_support = min_sample_support,
    support_fraction = support$support_fraction,
    n_supporting_samples = support$n_supporting_samples,
    n_total_samples = support$n_total_samples,
    recommended_global_cutoff = recommended_support$global_cutoff,
    recommended_support_fraction = recommended_support$support_fraction,
    recommended_n_supporting_samples = recommended_support$n_supporting_samples,
    min_cut = min_cut,
    max_cut = max_cut,
    bin_width = bin_width,
    min_drop_cells = min_drop_cells,
    diff_num = min_drop_cells,
    min_drop_fraction = min_drop_fraction,
    min_cells_after = min_cells_after,
    min_retention_after = min_retention_after,
    fallback_method = fallback_method,
    fallback_quantile = fallback_quantile
  )

  result <- list(
    cutoff = global_cutoff,
    sample_cutoff_summary = sample_cutoff_summary,
    retained_cell_profile = retained_cell_profile,
    interval_loss_table = interval_loss_table,
    significant_intervals = significant_intervals,
    settings = settings,
    sample_supported_global_cutoff = support,
    recommended_sample_supported_global_cutoff = recommended_support,
    counts_profile = retained_cell_profile_legacy,
    change_points = interval_loss_table
  )

  if (isTRUE(return_details)) {
    return(result)
  }

  global_cutoff
}

plot_scdetmito_outputs <- function(retained_cell_profile,
                                   interval_loss_table,
                                   sample_cutoff_summary,
                                   output_dir = ".") {
  output_dir <- ensure_output_dir(output_dir)
  sample_levels <- unique(retained_cell_profile$sample)
  palette <- scdetmito_palette(length(sample_levels))
  names(palette) <- sample_levels
  sample_cutoff_summary$sample_ordered <- stats::reorder(
    sample_cutoff_summary$sample,
    sample_cutoff_summary$detected_cutoff
  )
  cutoff_line_columns <- c(
    reference_cutoff = "Reference cutoff",
    first_significant_cutoff_high = "First significant boundary",
    largest_drop_cutoff = "Largest-drop cutoff",
    recommended_cutoff = "Recommended cutoff",
    selected_cutoff = "Selected cutoff"
  )
  cutoff_lines <- do.call(
    rbind,
    lapply(names(cutoff_line_columns), function(column) {
      if (!column %in% colnames(sample_cutoff_summary)) {
        return(NULL)
      }
      values <- as.numeric(sample_cutoff_summary[[column]])
      keep <- is.finite(values)
      if (!any(keep)) {
        return(NULL)
      }
      data.frame(
        sample = sample_cutoff_summary$sample[keep],
        cutoff = values[keep],
        cutoff_type = cutoff_line_columns[[column]],
        stringsAsFactors = FALSE
      )
    })
  )
  cutoff_lines_unique <- if (!is.null(cutoff_lines) && nrow(cutoff_lines)) {
    unique(cutoff_lines[, c("cutoff", "cutoff_type"), drop = FALSE])
  } else {
    data.frame(cutoff = numeric(), cutoff_type = character(), stringsAsFactors = FALSE)
  }
  cutoff_line_types <- c(
    "Reference cutoff" = "dashed",
    "First significant boundary" = "dashed",
    "Largest-drop cutoff" = "dashed",
    "Recommended cutoff" = "dotdash",
    "Selected cutoff" = "solid"
  )
  cutoff_line_colors <- c(
    "Reference cutoff" = "#6f7780",
    "First significant boundary" = "#2f6fb0",
    "Largest-drop cutoff" = "#b54b3a",
    "Recommended cutoff" = "#6a4ca3",
    "Selected cutoff" = "#1f1f1f"
  )
  add_cutoff_marker_layers <- function(plot_obj, marker_df, alpha_value = 0.7) {
    if (is.null(marker_df) || !nrow(marker_df)) {
      return(plot_obj)
    }
    marker_df$cutoff_type <- factor(marker_df$cutoff_type, levels = names(cutoff_line_types))
    for (marker_type in names(cutoff_line_types)) {
      type_df <- marker_df[marker_df$cutoff_type == marker_type, , drop = FALSE]
      if (!nrow(type_df)) {
        next
      }
      plot_obj <- plot_obj +
        ggplot2::geom_vline(
          data = type_df,
          ggplot2::aes(xintercept = .data$cutoff, linetype = .data$cutoff_type),
          inherit.aes = FALSE,
          color = cutoff_line_colors[[marker_type]],
          alpha = alpha_value,
          linewidth = if (marker_type %in% c("Selected cutoff", "Recommended cutoff")) 0.6 else 0.45,
          show.legend = TRUE
        )
    }
    plot_obj +
      ggplot2::scale_linetype_manual(
        values = cutoff_line_types,
        name = "Cutoff marker",
        breaks = names(cutoff_line_types)
      )
  }

  retained_plot <- ggplot2::ggplot(
    retained_cell_profile,
    ggplot2::aes(
      x = .data$cutoff,
      y = .data$retained_cells,
      color = .data$sample,
      group = .data$sample
    )
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::scale_x_reverse() +
    ggplot2::scale_color_manual(values = palette) +
    theme_scdetmito() +
    ggplot2::labs(
      title = "SCdetMito retained-cell profiles",
      x = "Mitochondrial ratio cutoff",
      y = "Retained cells",
      color = "Sample"
    )
  if (nrow(cutoff_lines_unique)) {
    retained_plot <- add_cutoff_marker_layers(retained_plot, cutoff_lines_unique, alpha_value = 0.55)
  }

  interval_plot <- ggplot2::ggplot(
    interval_loss_table,
    ggplot2::aes(
      x = .data$cutoff,
      y = .data$interval_loss,
      color = .data$is_significant,
      group = .data$sample
    )
  ) +
    ggplot2::geom_line(color = "#6d7882", linewidth = 0.45) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(~sample, scales = "free_y") +
    ggplot2::scale_x_reverse() +
    ggplot2::scale_color_manual(values = c("FALSE" = "#6d7882", "TRUE" = "#b54b3a")) +
    theme_scdetmito() +
    ggplot2::labs(
      title = "SCdetMito interval-specific retention-loss profiles",
      x = "Mitochondrial ratio cutoff",
      y = "Interval-specific cell loss",
      color = "Significant"
    )
  if (!is.null(cutoff_lines) && nrow(cutoff_lines)) {
    interval_plot <- add_cutoff_marker_layers(interval_plot, cutoff_lines, alpha_value = 0.72)
  }

  summary_plot <- ggplot2::ggplot(
    sample_cutoff_summary,
    ggplot2::aes(
      x = .data$sample_ordered,
      y = .data$detected_cutoff,
      fill = .data$cutoff_source,
      alpha = .data$cutoff_confidence
    )
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = scdetmito_palette(length(unique(sample_cutoff_summary$cutoff_source)))) +
    ggplot2::scale_alpha_manual(values = c("high" = 1, "moderate" = 0.78, "low" = 0.55, "user_defined" = 0.9)) +
    theme_scdetmito() +
    ggplot2::labs(
      title = "SCdetMito sample-level cutoff candidates",
      x = "Sample",
      y = "Detected mitochondrial ratio cutoff",
      fill = "Cutoff source",
      alpha = "Cutoff confidence"
    )

  ggplot2::ggsave(
    filename = file.path(output_dir, "SCdetMito_retained_cell_profiles.pdf"),
    plot = retained_plot,
    width = 10,
    height = 5,
    limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, "SCdetMito_interval_loss_profiles.pdf"),
    plot = interval_plot,
    width = 10,
    height = max(5, 2.6 * ceiling(length(sample_levels) / 2)),
    limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, "SCdetMito_sample_cutoff_summary.pdf"),
    plot = summary_plot,
    width = 8,
    height = max(4, 0.38 * nrow(sample_cutoff_summary) + 2),
    limitsize = FALSE
  )
}
