# SCdetMito
# Author: Silu Hu
# Contact: husilu0902@gmail.com
# Version: 1.3.6
# Last updated: 2026-05-20

# Internal helpers for QC workflows.

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1 && is.na(x))) {
    return(y)
  }
  x
}

ensure_output_dir <- function(output_dir = ".") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(output_dir, winslash = "/", mustWork = FALSE)
}

meta_column <- function(seurat_obj, column) {
  seurat_obj@meta.data[[column]]
}

theme_scdetmito <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#e7e7e7", linewidth = 0.25),
      strip.background = ggplot2::element_rect(fill = "#f2f3f3", color = "#c8c8c8"),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", color = "#1f2d35"),
      axis.title = ggplot2::element_text(color = "#24313a")
    )
}

scdetmito_palette <- function(n = 8) {
  base <- c(
    "#3b6f8f", "#b54b3a", "#4f8a5b", "#7b5ea7",
    "#c28f2c", "#4c8c8a", "#8d5b4c", "#68788f"
  )
  if (n <= length(base)) {
    return(base[seq_len(n)])
  }
  grDevices::colorRampPalette(base)(n)
}

validate_numeric_scalar <- function(value,
                                    name,
                                    lower = -Inf,
                                    upper = Inf,
                                    allow_infinite = FALSE) {
  if (!is.numeric(value) ||
    length(value) != 1 ||
    is.na(value) ||
    (!allow_infinite && !is.finite(value)) ||
    value < lower ||
    value > upper) {
    stop(
      "'", name, "' must be a numeric scalar in the interval [",
      lower, ", ", upper, "].",
      call. = FALSE
    )
  }

  value
}

normalize_mito_cutoff_value <- function(value,
                                        name = "max_mito",
                                        allow_scdet = FALSE) {
  if (allow_scdet && identical(value, "SCdetMito")) {
    return(value)
  }

  numeric_value <- suppressWarnings(as.numeric(value))
  if (length(numeric_value) != 1 || is.na(numeric_value) || !is.finite(numeric_value)) {
    expected <- if (isTRUE(allow_scdet)) {
      "a numeric mitochondrial fraction or 'SCdetMito'"
    } else {
      "a numeric mitochondrial fraction"
    }
    stop(
      "'", name, "' must be ", expected, ".",
      call. = FALSE
    )
  }

  if (numeric_value > 1 && numeric_value <= 100) {
    warning(
      "'", name, "' appears to be a percentage and was converted to a fraction.",
      call. = FALSE
    )
    numeric_value <- numeric_value / 100
  }

  if (numeric_value < 0 || numeric_value > 1) {
    stop(
      "'", name, "' must be between 0 and 1 when expressed as a fraction, ",
      "or between 0 and 100 when expressed as a percentage.",
      call. = FALSE
    )
  }

  numeric_value
}

validate_qc_bounds <- function(min_genes,
                               max_genes,
                               min_counts,
                               max_counts,
                               min_mito,
                               max_mito = NULL) {
  min_genes <- validate_numeric_scalar(min_genes, "min_genes", lower = 0)
  max_genes <- validate_numeric_scalar(max_genes, "max_genes", lower = 0, allow_infinite = TRUE)
  min_counts <- validate_numeric_scalar(min_counts, "min_counts", lower = 0)
  max_counts <- validate_numeric_scalar(max_counts, "max_counts", lower = 0, allow_infinite = TRUE)
  min_mito <- normalize_mito_cutoff_value(min_mito, "min_mito")
  max_mito <- if (is.null(max_mito)) {
    NULL
  } else {
    normalize_mito_cutoff_value(max_mito, "max_mito")
  }

  if (min_genes > max_genes) {
    stop("'min_genes' cannot be greater than 'max_genes'.", call. = FALSE)
  }
  if (min_counts > max_counts) {
    stop("'min_counts' cannot be greater than 'max_counts'.", call. = FALSE)
  }
  if (!is.null(max_mito) && min_mito > max_mito) {
    stop("'min_mito' cannot be greater than 'max_mito'.", call. = FALSE)
  }

  list(
    min_genes = min_genes,
    max_genes = max_genes,
    min_counts = min_counts,
    max_counts = max_counts,
    min_mito = min_mito,
    max_mito = max_mito
  )
}

build_qc_summary <- function(seurat_obj,
                             nFeature_RNA = "nFeature_RNA",
                             nCount_RNA = "nCount_RNA",
                             mitoRatio = "mitoRatio") {
  data.frame(
    CellCounts = ncol(seurat_obj),
    GenesCounts_median = stats::median(meta_column(seurat_obj, nFeature_RNA), na.rm = TRUE),
    TranscriptCounts_median = stats::median(meta_column(seurat_obj, nCount_RNA), na.rm = TRUE),
    MitoRatio_median = stats::median(meta_column(seurat_obj, mitoRatio), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

add_qc_outcome_metrics <- function(summary_df, reference_cell_count) {
  summary_df$CellsRemoved <- reference_cell_count - summary_df$CellCounts
  summary_df$RetentionRate <- round(summary_df$CellCounts / reference_cell_count, 4)
  summary_df
}

filter_cells_by_metrics <- function(seurat_obj,
                                    nFeature_RNA = "nFeature_RNA",
                                    nCount_RNA = "nCount_RNA",
                                    mitoRatio = "mitoRatio",
                                    min_genes = 200,
                                    max_genes = Inf,
                                    min_counts = 500,
                                    max_counts = Inf,
                                    min_mito = 0,
                                    max_mito = 1,
                                    max_mito_by_cell = NULL) {
  qc_bounds <- validate_qc_bounds(
    min_genes = min_genes,
    max_genes = max_genes,
    min_counts = min_counts,
    max_counts = max_counts,
    min_mito = min_mito,
    max_mito = max_mito
  )
  min_genes <- qc_bounds$min_genes
  max_genes <- qc_bounds$max_genes
  min_counts <- qc_bounds$min_counts
  max_counts <- qc_bounds$max_counts
  min_mito <- qc_bounds$min_mito
  max_mito <- qc_bounds$max_mito

  metadata <- seurat_obj@meta.data
  if (is.null(max_mito_by_cell)) {
    mito_upper <- rep(max_mito, nrow(metadata))
  } else if (is.null(names(max_mito_by_cell))) {
    if (length(max_mito_by_cell) != nrow(metadata)) {
      stop("'max_mito_by_cell' must have one entry per cell.", call. = FALSE)
    }
    mito_upper <- as.numeric(max_mito_by_cell)
  } else {
    mito_upper <- as.numeric(max_mito_by_cell[rownames(metadata)])
  }

  if (anyNA(mito_upper)) {
    stop("Cell-level mitochondrial cutoffs could not be resolved for all cells.", call. = FALSE)
  }

  keep <- metadata[[nFeature_RNA]] >= min_genes &
    metadata[[nFeature_RNA]] <= max_genes &
    metadata[[nCount_RNA]] >= min_counts &
    metadata[[nCount_RNA]] <= max_counts &
    metadata[[mitoRatio]] >= min_mito &
    metadata[[mitoRatio]] <= mito_upper

  keep[is.na(keep)] <- FALSE
  kept_cells <- rownames(metadata)[keep]

  if (!length(kept_cells)) {
    stop("No cells passed the requested QC thresholds.", call. = FALSE)
  }

  subset(seurat_obj, cells = kept_cells)
}

run_loss_enrichment_test <- function(loss_value,
                                     background_losses,
                                     method = c("mad_zscore", "empirical_tail", "poisson_tail", "zscore", "threshold_only")) {
  method <- match.arg(method)
  background_losses <- as.numeric(background_losses)
  background_losses <- background_losses[is.finite(background_losses) & background_losses >= 0]

  if (identical(method, "threshold_only")) {
    return(list(
      p_value = NA_real_,
      loss_score = as.numeric(loss_value),
      test_method = "threshold_only",
      baseline_loss = if (length(background_losses)) mean(background_losses, na.rm = TRUE) else NA_real_
    ))
  }

  if (!length(background_losses)) {
    return(list(
      p_value = 1,
      loss_score = NA_real_,
      test_method = paste0(method, "_insufficient_baseline"),
      baseline_loss = NA_real_
    ))
  }

  baseline_loss <- mean(background_losses, na.rm = TRUE)
  if (!is.finite(baseline_loss)) {
    return(list(
      p_value = 1,
      loss_score = NA_real_,
      test_method = paste0(method, "_insufficient_baseline"),
      baseline_loss = NA_real_
    ))
  }

  if (loss_value <= 0) {
    return(list(
      p_value = 1,
      loss_score = 0,
      test_method = method,
      baseline_loss = if (identical(method, "mad_zscore")) {
        stats::median(background_losses, na.rm = TRUE)
      } else {
        baseline_loss
      }
    ))
  }

  robust_loss_reference <- function(values) {
    center <- stats::median(values, na.rm = TRUE)
    spread <- stats::mad(values, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(spread) || spread <= 0) {
      iqr_spread <- stats::IQR(values, na.rm = TRUE) / 1.349
      if (is.finite(iqr_spread) && iqr_spread > 0) {
        spread <- iqr_spread
      }
    }
    if (!is.finite(spread) || spread <= 0) {
      sd_spread <- stats::sd(values, na.rm = TRUE)
      if (is.finite(sd_spread) && sd_spread > 0) {
        spread <- sd_spread
      }
    }
    if (!is.finite(spread) || spread <= 0) {
      spread <- .Machine$double.eps
    }
    list(center = center, spread = spread)
  }

  loss_score <- switch(method,
    empirical_tail = loss_value,
    poisson_tail = loss_value / max(baseline_loss, 1e-8),
    zscore = {
      background_sd <- stats::sd(background_losses, na.rm = TRUE)
      if (!is.finite(background_sd) || background_sd <= 0) {
        if (loss_value > baseline_loss) Inf else 0
      } else {
        (loss_value - baseline_loss) / background_sd
      }
    },
    mad_zscore = {
      robust_ref <- robust_loss_reference(background_losses)
      (loss_value - robust_ref$center) / robust_ref$spread
    }
  )

  p_value <- switch(method,
    empirical_tail = {
      (sum(background_losses >= loss_value, na.rm = TRUE) + 1) /
        (length(background_losses) + 1)
    },
    poisson_tail = stats::ppois(
      q = loss_value - 1,
      lambda = max(baseline_loss, 1e-8),
      lower.tail = FALSE
    ),
    zscore = {
      background_sd <- stats::sd(background_losses, na.rm = TRUE)
      if (!is.finite(background_sd) || background_sd <= 0) {
        if (loss_value > baseline_loss) .Machine$double.xmin else 1
      } else {
        z_score <- (loss_value - baseline_loss) / background_sd
        stats::pnorm(z_score, lower.tail = FALSE)
      }
    },
    mad_zscore = {
      robust_ref <- robust_loss_reference(background_losses)
      z_score <- (loss_value - robust_ref$center) / robust_ref$spread
      stats::pnorm(z_score, lower.tail = FALSE)
    }
  )

  list(
    p_value = max(p_value, .Machine$double.xmin),
    loss_score = loss_score,
    test_method = method,
    baseline_loss = if (identical(method, "mad_zscore")) {
      stats::median(background_losses, na.rm = TRUE)
    } else {
      baseline_loss
    }
  )
}

select_detected_cutoff <- function(group_df,
                                   fallback_cutoff,
                                   alpha = 0.05,
                                   method = c("largest_drop", "max_significant", "median_significant", "min_significant"),
                                   snap_fun = identity) {
  method <- match.arg(method)
  if (!nrow(group_df)) {
    return(list(
      cutoff = fallback_cutoff,
      source = "no_candidates",
      significant_cutoff_count = 0L
    ))
  }

  significant_col <- if ("is_significant" %in% colnames(group_df)) {
    "is_significant"
  } else {
    NULL
  }
  significant_df <- if (!is.null(significant_col)) {
    group_df[group_df[[significant_col]] %in% TRUE, , drop = FALSE]
  } else {
    group_df[group_df$p_adj <= alpha, , drop = FALSE]
  }
  if (!nrow(significant_df)) {
    return(list(
      cutoff = fallback_cutoff,
      source = "no_significant_candidates",
      significant_cutoff_count = 0L
    ))
  }

  selected_cutoff <- switch(method,
    largest_drop = {
      ranked_df <- significant_df[order(-significant_df$interval_loss, -significant_df$cutoff), , drop = FALSE]
      ranked_df$cutoff[[1]]
    },
    max_significant = max(significant_df$cutoff, na.rm = TRUE),
    median_significant = {
      snap_fun(stats::median(unique(significant_df$cutoff), na.rm = TRUE))
    },
    min_significant = min(significant_df$cutoff, na.rm = TRUE)
  )

  list(
    cutoff = snap_fun(selected_cutoff),
    source = paste0("significant_", method),
    significant_cutoff_count = length(unique(significant_df$cutoff))
  )
}

classify_cutoff_confidence <- function(fallback_used,
                                       cutoff_source,
                                       significant_interval_count,
                                       largest_interval_loss_fraction) {
  if (isTRUE(fallback_used) || identical(cutoff_source, "not_detected")) {
    return("low")
  }
  if (!is.finite(largest_interval_loss_fraction)) {
    return("low")
  }
  if (significant_interval_count >= 2L && largest_interval_loss_fraction >= 0.02) {
    return("high")
  }
  if (significant_interval_count >= 1L && largest_interval_loss_fraction >= 0.01) {
    return("moderate")
  }
  "low"
}

classify_detector_mode <- function(loss_test, p_adj_method) {
  if (identical(loss_test, "threshold_only")) {
    return("rule_based_sensitivity")
  }
  if (identical(loss_test, "mad_zscore") && !identical(p_adj_method, "none")) {
    return("primary_fdr_mode")
  }
  if (identical(p_adj_method, "none")) {
    return("unadjusted_sensitivity_mode")
  }
  "sensitivity_exploratory_mode"
}

capture_session_info <- function() {
  paste(utils::capture.output(utils::sessionInfo()), collapse = "\n")
}

package_version_or_na <- function(pkg = "SCdetMito") {
  version <- tryCatch(
    as.character(utils::packageVersion(pkg)),
    error = function(e) NA_character_
  )
  version
}

build_qc_provenance <- function(function_name,
                                input_cell_count,
                                output_cell_count,
                                parameters = list(),
                                cutoff_plan = NULL,
                                cutoff_source = NULL,
                                cutoff_confidence = NULL,
                                loss_test = NULL,
                                sample_cutoff_method = NULL,
                                warnings = character(),
                                fallback_events = NULL,
                                doublet_filtering = NULL) {
  output_cell_count <- as.numeric(output_cell_count)
  input_cell_count <- as.numeric(input_cell_count)
  list(
    package = "SCdetMito",
    package_version = package_version_or_na("SCdetMito"),
    function_name = function_name,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    input_cell_count = input_cell_count,
    output_cell_count = output_cell_count,
    retention_rate = if (is.finite(input_cell_count) && input_cell_count > 0) {
      output_cell_count / input_cell_count
    } else {
      NA_real_
    },
    parameters = parameters,
    cutoff_plan = cutoff_plan,
    cutoff_source = cutoff_source,
    cutoff_confidence = cutoff_confidence,
    loss_test = loss_test,
    sample_cutoff_method = sample_cutoff_method,
    sessionInfo = capture_session_info(),
    warnings = warnings,
    fallback_events = fallback_events,
    doublet_filtering = doublet_filtering %||% list(
      requested = FALSE,
      status = "skipped",
      error_message = NA_character_
    )
  )
}

attach_qc_provenance_to_cutoff_plan <- function(cutoff_plan, provenance) {
  qc_misc <- cutoff_plan
  for (field in setdiff(names(provenance), "cutoff_plan")) {
    qc_misc[[field]] <- provenance[[field]]
  }
  qc_misc$provenance <- provenance
  qc_misc
}

run_doublet_filter_safely <- function(seurat_obj,
                                      requested = FALSE,
                                      fail_action = c("warn", "stop")) {
  fail_action <- match.arg(fail_action)
  status <- list(
    requested = isTRUE(requested),
    status = if (isTRUE(requested)) "skipped" else "not_requested",
    error_message = NA_character_
  )
  if (!isTRUE(requested)) {
    return(list(object = seurat_obj, status = status))
  }

  result <- tryCatch(
    {
      filtered <- run_doublet_filter(seurat_obj)
      status$status <- "completed"
      list(object = filtered, status = status)
    },
    error = function(e) {
      status$status <- "failed"
      status$error_message <- conditionMessage(e)
      if (identical(fail_action, "stop")) {
        stop(e)
      }
      warning(
        "DoubletFinder filtering failed or is unavailable; returning the mito-QC-filtered object. ",
        "Reason: ", conditionMessage(e),
        call. = FALSE
      )
      list(object = seurat_obj, status = status)
    },
    warning = function(w) {
      status$status <- "skipped"
      status$error_message <- conditionMessage(w)
      if (identical(fail_action, "stop")) {
        stop(w)
      }
      warning(
        "DoubletFinder filtering was skipped; returning the mito-QC-filtered object. ",
        "Reason: ", conditionMessage(w),
        call. = FALSE
      )
      list(object = seurat_obj, status = status)
    }
  )
  result
}

run_doublet_filter <- function(seurat_obj) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    stop(
      "Doublet removal requires the optional package 'DoubletFinder'. ",
      "Install it separately or call the QC function with removeDouble = FALSE.",
      call. = FALSE
    )
  }

  if (ncol(seurat_obj) < 50) {
    warning(
      "Doublet removal was skipped because fewer than 50 cells remained after QC filtering.",
      call. = FALSE
    )
    return(seurat_obj)
  }

  seu <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(
    object = seu,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  seu <- Seurat::ScaleData(object = seu, verbose = FALSE)
  max_pc <- min(30L, max(5L, ncol(seu) - 1L, na.rm = TRUE))
  seu <- Seurat::RunPCA(
    object = seu,
    features = Seurat::VariableFeatures(object = seu),
    npcs = max_pc,
    verbose = FALSE
  )
  pc.num <- seq_len(max_pc)

  sweep.res.list <- DoubletFinder::paramSweep(
    seu,
    PCs = pc.num,
    sct = FALSE
  )
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  pK_bcmvn <- as.numeric(as.character(pK_bcmvn))

  if (!length(pK_bcmvn) || is.na(pK_bcmvn)) {
    warning("DoubletFinder could not determine a preferred pK; skipping doublet removal.", call. = FALSE)
    return(seurat_obj)
  }

  doublet_rate <- ncol(seu) / 5000 * 0.039
  nExp_poi <- round(doublet_rate * ncol(seu))

  if (nExp_poi < 1) {
    warning("Estimated expected doublets was below 1; skipping doublet removal.", call. = FALSE)
    return(seurat_obj)
  }

  seu <- DoubletFinder::doubletFinder(
    seu,
    PCs = pc.num,
    pN = 0.25,
    pK = pK_bcmvn,
    nExp = nExp_poi,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  colnames(seu@meta.data)[ncol(seu@meta.data)] <- "double_info"
  singlet_cells <- rownames(seu@meta.data)[seu@meta.data$double_info == "Singlet"]
  subset(seu, cells = singlet_cells)
}
