# SCdetMito
# Author: Silu Hu
# Contact: husilu0902@gmail.com
# Version: 1.4.3
# Last updated: 2026-05-23

# Avoid R CMD check notes for ggplot2 non-standard evaluation columns.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "parameter_label",
    "sample",
    "detected_cutoff",
    "loss_test",
    "cutoff_median",
    "fallback_fraction",
    "not_detected_fraction",
    "cutoff_min",
    "cutoff_max"
  ))
}

#' Sensitivity analysis for SCdetMito cutoff estimation
#'
#' @description
#' `SCdetMito_sensitivity()` evaluates whether sample-level mitochondrial
#' cutoffs are stable across user-specified retention-loss enrichment settings.
#' It is intended for method auditing, supplementary analyses, and manuscript
#' robustness checks.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq QC metadata.
#' @param sample_col,by Metadata column defining sample IDs. `sample_col` is the
#'   preferred name; `by` is retained for backward compatibility.
#' @param mito_col,mitoRatio Metadata column storing mitochondrial ratios.
#' @param bin_width Candidate cutoff grid step sizes to test.
#' @param alpha Adjusted p-value cutoffs to test.
#' @param min_drop_fraction Minimum interval-loss fractions to test.
#' @param loss_test Loss enrichment tests to evaluate. Defaults to
#'   `"mad_zscore"` plus `"empirical_tail"` as a sensitivity comparison.
#'   Supported values are `"mad_zscore"`, `"empirical_tail"`,
#'   `"poisson_tail"`, `"zscore"`, and `"threshold_only"`.
#' @param sample_cutoff_method Sample-level cutoff selection rules to evaluate.
#' @param write_tables Whether to write the combined sensitivity table.
#' @param write_plots Whether to write optional stability plots.
#' @param output_dir Output directory for optional tables and plots.
#' @param return_details Whether to return both the per-combination table and
#'   sample-level sensitivity summary. Defaults to `FALSE`, which preserves the
#'   legacy behavior of returning only the per-combination table.
#' @param ... Additional arguments forwarded to [SCdetMito()].
#'
#' @return A data frame with one row per sample and parameter combination, or a
#'   list with `sensitivity_table` and `sensitivity_summary` when
#'   `return_details = TRUE`.
#' @export
SCdetMito_sensitivity <- function(seurat_obj,
                                  sample_col = NULL,
                                  by = "sample",
                                  mito_col = NULL,
                                  mitoRatio = "mitoRatio",
                                  bin_width = c(0.01, 0.02),
                                  alpha = c(0.01, 0.05, 0.10),
                                  min_drop_fraction = c(0.005, 0.01, 0.02),
                                  loss_test = c("mad_zscore", "empirical_tail"),
                                  sample_cutoff_method = c("largest_drop", "first_significant_high", "reference_guided"),
                                  write_tables = FALSE,
                                  write_plots = FALSE,
                                  output_dir = ".",
                                  return_details = FALSE,
                                  ...) {
  if (!is.null(sample_col)) {
    by <- sample_col
  }
  if (!is.null(mito_col)) {
    mitoRatio <- mito_col
  }
  if (is.null(by) || !nzchar(by)) {
    stop("'sample_col' or backward-compatible 'by' must be provided.", call. = FALSE)
  }

  loss_test <- match.arg(
    loss_test,
    choices = c("mad_zscore", "empirical_tail", "poisson_tail", "zscore", "threshold_only"),
    several.ok = TRUE
  )
  sample_cutoff_method <- match.arg(
    sample_cutoff_method,
    choices = c(
      "largest_drop",
      "first_significant_high",
      "first_significant_low",
      "reference_guided",
      "largest_drop_with_reference_guard",
      "max_significant",
      "median_significant",
      "min_significant"
    ),
    several.ok = TRUE
  )

  grid <- expand.grid(
    bin_width = bin_width,
    alpha = alpha,
    min_drop_fraction = min_drop_fraction,
    loss_test = loss_test,
    sample_cutoff_method = sample_cutoff_method,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  rows <- vector("list", nrow(grid))
  for (idx in seq_len(nrow(grid))) {
    combo <- grid[idx, , drop = FALSE]
    parameter_id <- paste0("P", sprintf("%03d", idx))
    detection <- SCdetMito(
      seurat_obj = seurat_obj,
      by = by,
      mitoRatio = mitoRatio,
      bin_width = combo$bin_width,
      alpha = combo$alpha,
      min_drop_fraction = combo$min_drop_fraction,
      loss_test = combo$loss_test,
      sample_cutoff_method = combo$sample_cutoff_method,
      write_tables = FALSE,
      write_plots = FALSE,
      plot = FALSE,
      table_out = FALSE,
      return_details = TRUE,
      ...
    )

    sample_summary <- detection$sample_cutoff_summary
    sample_summary$parameter_id <- parameter_id
    sample_summary$bin_width <- combo$bin_width
    sample_summary$alpha <- combo$alpha
    sample_summary$min_drop_fraction <- combo$min_drop_fraction
    sample_summary$loss_test <- combo$loss_test
    sample_summary$sample_cutoff_method <- combo$sample_cutoff_method
    sample_summary$global_cutoff <- detection$cutoff
    rows[[idx]] <- sample_summary
  }

  sensitivity_table <- do.call(rbind, rows)
  requested_cols <- c(
    "sample",
    "parameter_id",
    "bin_width",
    "alpha",
    "min_drop_fraction",
    "loss_test",
    "sample_cutoff_method",
    "detected_cutoff",
    "cutoff_source",
    "retained_cells_at_cutoff",
    "retention_fraction_at_cutoff",
    "global_cutoff"
  )
  sensitivity_table <- sensitivity_table[, c(
    requested_cols,
    setdiff(colnames(sensitivity_table), requested_cols)
  ), drop = FALSE]
  sensitivity_summary <- summarize_scdetmito_sensitivity(sensitivity_table)

  if (isTRUE(write_tables) || isTRUE(write_plots)) {
    output_dir <- ensure_output_dir(output_dir)
  }
  if (isTRUE(write_tables)) {
    utils::write.csv(
      sensitivity_table,
      file = file.path(output_dir, "SCdetMito_sensitivity_results.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      sensitivity_summary,
      file = file.path(output_dir, "SCdetMito_sensitivity_summary.csv"),
      row.names = FALSE
    )
  }
  if (isTRUE(write_plots)) {
    write_scdetmito_sensitivity_plots(
      sensitivity_table = sensitivity_table,
      sensitivity_summary = sensitivity_summary,
      output_dir = output_dir
    )
  }

  if (isTRUE(return_details)) {
    return(list(
      sensitivity_table = sensitivity_table,
      sensitivity_summary = sensitivity_summary
    ))
  }
  sensitivity_table
}

summarize_scdetmito_sensitivity <- function(sensitivity_table) {
  sample_levels <- unique(as.character(sensitivity_table$sample))
  rows <- lapply(sample_levels, function(sample_id) {
    sample_df <- sensitivity_table[as.character(sensitivity_table$sample) == sample_id, , drop = FALSE]
    finite_cutoffs <- sample_df$detected_cutoff[is.finite(sample_df$detected_cutoff)]
    quantiles <- if (length(finite_cutoffs)) {
      stats::quantile(finite_cutoffs, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE, type = 8)
    } else {
      c(NA_real_, NA_real_, NA_real_)
    }
    data.frame(
      sample = sample_id,
      n_parameter_combinations = nrow(sample_df),
      cutoff_median = quantiles[[2]],
      cutoff_iqr = quantiles[[3]] - quantiles[[1]],
      cutoff_min = if (length(finite_cutoffs)) min(finite_cutoffs, na.rm = TRUE) else NA_real_,
      cutoff_max = if (length(finite_cutoffs)) max(finite_cutoffs, na.rm = TRUE) else NA_real_,
      fallback_fraction = mean(sample_df$fallback_used %in% TRUE, na.rm = TRUE),
      not_detected_fraction = mean(sample_df$cutoff_source == "not_detected", na.rm = TRUE),
      low_confidence_fraction = if ("cutoff_confidence" %in% colnames(sample_df)) {
        mean(sample_df$cutoff_confidence == "low", na.rm = TRUE)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

write_scdetmito_sensitivity_plots <- function(sensitivity_table,
                                              sensitivity_summary = NULL,
                                              output_dir = ".") {
  plot_df <- sensitivity_table
  plot_df$parameter_label <- paste(
    plot_df$loss_test,
    plot_df$sample_cutoff_method,
    paste0("bw=", plot_df$bin_width),
    paste0("alpha=", plot_df$alpha),
    paste0("drop=", plot_df$min_drop_fraction),
    sep = "\n"
  )

  heatmap_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = parameter_label, y = sample, fill = detected_cutoff)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(low = "#f6f7fb", high = "#274472", na.value = "#d9dde7") +
    theme_scdetmito(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    ggplot2::labs(x = "Parameter combination", y = "Sample", fill = "Detected cutoff")

  distribution_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = detected_cutoff, fill = loss_test)
  ) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.75, position = "identity") +
    ggplot2::scale_fill_manual(values = scdetmito_palette(length(unique(plot_df$loss_test)))) +
    theme_scdetmito() +
    ggplot2::labs(x = "Detected sample-level cutoff", y = "Count", fill = "Loss test")

  ggplot2::ggsave(
    filename = file.path(output_dir, "SCdetMito_sensitivity_cutoff_stability_heatmap.pdf"),
    plot = heatmap_plot,
    width = 10,
    height = max(4, min(12, length(unique(plot_df$sample)) * 0.35 + 2))
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, "SCdetMito_sensitivity_cutoff_distribution.pdf"),
    plot = distribution_plot,
    width = 7,
    height = 4.8
  )

  if (!is.null(sensitivity_summary) && nrow(sensitivity_summary)) {
    summary_plot <- ggplot2::ggplot(
      sensitivity_summary,
      ggplot2::aes(x = stats::reorder(sample, cutoff_median), y = cutoff_median)
    ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = cutoff_min, ymax = cutoff_max),
        width = 0.18,
        color = "#5b6b7b"
      ) +
      ggplot2::geom_point(ggplot2::aes(color = fallback_fraction), size = 2.5) +
      ggplot2::coord_flip() +
      ggplot2::scale_color_gradient(low = "#2f7d4a", high = "#b33a3a") +
      theme_scdetmito() +
      ggplot2::labs(
        x = "Sample",
        y = "Median detected cutoff across sensitivity grid",
        color = "Fallback fraction"
      )
    ggplot2::ggsave(
      filename = file.path(output_dir, "SCdetMito_sensitivity_summary.pdf"),
      plot = summary_plot,
      width = 7,
      height = max(4, min(12, nrow(sensitivity_summary) * 0.35 + 2))
    )
  }
}
