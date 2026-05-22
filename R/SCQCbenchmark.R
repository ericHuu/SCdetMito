# SCdetMito
# Author: Silu Hu
# Contact: husilu0902@gmail.com
# Version: 1.4.1
# Last updated: 2026-05-21

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "strategy_order",
    "total_cells_after",
    "cutoff_strategy",
    "median_mito_after",
    "level_id",
    "strategy",
    "retention_rate",
    "overall_score",
    "scoring_profile"
  ))
}

#' SCQCbenchmark: benchmark fixed and adaptive QC strategies
#'
#' @description
#' `SCQCbenchmark()` compares fixed mitochondrial thresholds with adaptive
#' `SCdetMito`-based strategies on the same multi-sample Seurat object.
#'
#' @details
#' The function is intended for method evaluation and manuscript preparation.
#' It summarizes retained cell numbers, residual mitochondrial burden,
#' sample/group balance, and lightweight downstream clustering diagnostics for
#' each QC strategy. Composite scores are user-adjustable decision-support
#' summaries and do not imply a universal QC optimum.
#'
#' @param seurat_obj A Seurat object containing multi-sample single-cell data.
#' @param sample_by,sample_col Metadata column defining sample IDs.
#'   `sample_col` is the clearer preferred name; `sample_by` is retained for
#'   backward compatibility.
#' @param group_by,group_col Optional metadata column defining biological groups.
#'   `group_col` is the clearer preferred name; `group_by` is retained for
#'   backward compatibility.
#' @param strategies Optional list of strategy specifications. When `NULL`, a
#'   default benchmark panel is used.
#' @param scdet_options Optional list of arguments forwarded to [SCdetMito()]
#'   for adaptive strategies. Individual strategies can override these values
#'   by defining their own `scdet_options` entry.
#' @param nFeature_RNA,feature_col Metadata column storing detected feature
#'   counts.
#' @param nCount_RNA,count_col Metadata column storing UMI/count totals.
#' @param mitoRatio,mito_col Metadata column storing mitochondrial ratios.
#' @param min_genes Minimum number of detected genes. Defaults to `200`.
#' @param max_genes Maximum number of detected genes. Defaults to `Inf`.
#' @param min_counts Minimum number of counts. Defaults to `500`.
#' @param max_counts Maximum number of counts. Defaults to `Inf`.
#' @param min_mito Minimum mitochondrial ratio. Defaults to `0`. Values between
#'   `1` and `100` are interpreted as percentages and converted to fractions by
#'   the underlying QC functions.
#' @param removeDouble,remove_doublets Whether to run optional DoubletFinder-
#'   based doublet removal inside each benchmark strategy. Defaults to `FALSE`.
#' @param run_downstream Whether to run lightweight downstream clustering and
#'   compute cluster-mixing diagnostics. Defaults to `TRUE`.
#' @param cluster_resolution Resolution used for Seurat clustering in downstream
#'   benchmarking. Defaults to `0.5`.
#' @param npcs Maximum number of PCs used for downstream clustering. Defaults to
#'   `20`.
#' @param seed Random seed for downstream analysis. Defaults to `123`.
#' @param score_weights Optional named numeric vector or list controlling a
#'   custom composite score. When supplied, it is reported as the `custom`
#'   scoring profile and used as the backward-compatible primary score.
#' @param scoring_profiles Built-in scoring profiles to compute. Available
#'   profiles are `"balanced"`, `"retention_focused"`, `"stringent_mito"`, and
#'   `"balance_focused"`.
#' @param species,tissue Optional species and tissue names forwarded to
#'   [SCdetMito()] for reference-aware adaptive strategies.
#' @param reference_cutoff,reference_table,reference_warning Optional reference
#'   cutoff controls forwarded to [SCdetMito()].
#' @param plot,write_plots Whether to export benchmark plots. Defaults to
#'   `TRUE`. `write_plots` is the clearer preferred name.
#' @param save_objects Whether to save strategy-specific filtered Seurat objects.
#'   Defaults to `FALSE`.
#' @param output_dir Directory for benchmark outputs. Defaults to `"."`.
#' @param auto_add_mito Whether to calculate `mito_col` with
#'   [ensure_mito_ratio()] when the column is absent. Defaults to `TRUE`.
#' @param mito_features Optional complete mitochondrial feature vector forwarded
#'   to [ensure_mito_ratio()] and adaptive SCdetMito strategies.
#' @param mito_pattern Optional mitochondrial feature regex forwarded to
#'   [ensure_mito_ratio()] and adaptive SCdetMito strategies.
#' @param mito_assay Optional assay used for mitochondrial ratio calculation.
#' @param recompute_mito Whether to recompute `mito_col` even when it already
#'   exists. Defaults to `FALSE`.
#'
#' @return Invisibly returns a list with raw benchmark tables, normalized score
#'   tables, weights used, strategies recommended under each scoring profile,
#'   retention tables, sample-level cutoff plans, and filtered objects. Legacy
#'   aliases `summary`, `scores`, and `recommendation` are retained.
#' @export
SCQCbenchmark <- function(seurat_obj,
                          sample_by = NULL,
                          group_by = NULL,
                          sample_col = NULL,
                          group_col = NULL,
                          strategies = NULL,
                          scdet_options = list(),
                          nFeature_RNA = "nFeature_RNA",
                          nCount_RNA = "nCount_RNA",
                          mitoRatio = "mitoRatio",
                          feature_col = NULL,
                          count_col = NULL,
                          mito_col = NULL,
                          min_genes = 200,
                          max_genes = Inf,
                          min_counts = 500,
                          max_counts = Inf,
                          min_mito = 0,
                          remove_doublets = NULL,
                          removeDouble = FALSE,
                          run_downstream = TRUE,
                          cluster_resolution = 0.5,
                          npcs = 20,
                          seed = 123,
                          score_weights = NULL,
                          scoring_profiles = c(
                            "balanced",
                            "retention_focused",
                            "stringent_mito",
                            "balance_focused"
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
                          write_plots = NULL,
                          plot = TRUE,
                          save_objects = FALSE,
                          output_dir = ".") {
  output_dir <- ensure_output_dir(output_dir)
  if (!is.null(sample_col)) {
    sample_by <- sample_col
  }
  if (is.null(sample_by) || !nzchar(sample_by)) {
    stop("'sample_col' or backward-compatible 'sample_by' must be provided.", call. = FALSE)
  }
  if (!is.null(group_col)) {
    group_by <- group_col
  }
  group_by <- normalize_group_by(group_by)
  if (!is.null(feature_col)) {
    nFeature_RNA <- feature_col
  }
  if (!is.null(count_col)) {
    nCount_RNA <- count_col
  }
  if (!is.null(mito_col)) {
    mitoRatio <- mito_col
  }
  if (!is.null(remove_doublets)) {
    removeDouble <- remove_doublets
  }
  if (!is.null(write_plots)) {
    plot <- write_plots
  }

  seurat_obj <- check_seu(seurat_obj, sample_by)
  if (!is.null(group_by)) {
    seurat_obj <- check_seu(seurat_obj, group_by)
  }
  seurat_obj <- check_seu(seurat_obj, nFeature_RNA, must_be_numeric = TRUE)
  seurat_obj <- check_seu(seurat_obj, nCount_RNA, must_be_numeric = TRUE)
  if (!mitoRatio %in% colnames(seurat_obj@meta.data) && !isTRUE(auto_add_mito)) {
    stop(
      "Mitochondrial column '",
      mitoRatio,
      "' is missing and auto_add_mito = FALSE. Provide a valid mito_col or set auto_add_mito = TRUE.",
      call. = FALSE
    )
  }
  seurat_obj <- ensure_mito_ratio(
    object = seurat_obj,
    mito_col = mitoRatio,
    species = species,
    assay = mito_assay,
    mito_features = mito_features,
    mito_pattern = mito_pattern,
    overwrite = isTRUE(recompute_mito),
    convert_percent = TRUE,
    verbose = FALSE
  )
  if (!is.list(scdet_options)) {
    stop("'scdet_options' must be a list.", call. = FALSE)
  }
  scoring_profile_list <- resolve_scoring_profiles(
    scoring_profiles = scoring_profiles,
    score_weights = score_weights
  )
  scdet_options <- utils::modifyList(
    list(
      species = species,
      tissue = tissue,
      reference_cutoff = reference_cutoff,
      reference_table = reference_table,
      reference_warning = reference_warning,
      auto_add_mito = auto_add_mito,
      mito_features = mito_features,
      mito_pattern = mito_pattern,
      mito_assay = mito_assay,
      recompute_mito = recompute_mito
    ),
    scdet_options
  )
  primary_profile <- names(scoring_profile_list)[1]

  strategies <- normalize_benchmark_strategies(strategies, group_by = group_by)
  strategy_names <- vapply(strategies, `[[`, character(1), "name")
  if (anyDuplicated(strategy_names)) {
    stop("Strategy names must be unique.", call. = FALSE)
  }

  summary_rows <- list()
  sample_retention_rows <- list()
  group_retention_rows <- list()
  cutoff_plan_rows <- list()
  filtered_objects <- list()

  for (strategy in strategies) {
    strategy_dir <- ensure_output_dir(file.path(output_dir, strategy$name))
    strategy_scdet_options <- utils::modifyList(
      scdet_options,
      strategy$scdet_options
    )
    qc_obj <- suppressWarnings(SCQCmulti(
      seurat_obj = seurat_obj,
      by = sample_by,
      group_by = group_by,
      mode = strategy$mode,
      cutoff_strategy = strategy$cutoff_strategy,
      scdet_options = strategy_scdet_options,
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      mitoRatio = mitoRatio,
      min_genes = min_genes,
      max_genes = max_genes,
      min_counts = min_counts,
      max_counts = max_counts,
      min_mito = min_mito,
      max_mito = strategy$max_mito,
      removeDouble = removeDouble,
      plot = FALSE,
      table_out = TRUE,
      output_dir = strategy_dir
    ))

    summary_rows[[strategy$name]] <- summarize_benchmark_strategy(
      original_obj = seurat_obj,
      filtered_obj = qc_obj,
      strategy = strategy,
      sample_by = sample_by,
      group_by = group_by,
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      mitoRatio = mitoRatio,
      run_downstream = run_downstream,
      cluster_resolution = cluster_resolution,
      npcs = npcs,
      seed = seed
    )
    sample_retention_rows[[strategy$name]] <- build_retention_table(
      original_obj = seurat_obj,
      filtered_obj = qc_obj,
      column = sample_by,
      level_type = "sample",
      strategy_name = strategy$name
    )

    if (!is.null(group_by)) {
      group_retention_rows[[strategy$name]] <- build_retention_table(
        original_obj = seurat_obj,
        filtered_obj = qc_obj,
        column = group_by,
        level_type = "group",
        strategy_name = strategy$name
      )
    }

    if (!is.null(qc_obj@misc$SCdetMito_QC$sample_plan)) {
      plan_df <- qc_obj@misc$SCdetMito_QC$sample_plan
      plan_df$strategy <- strategy$name
      cutoff_plan_rows[[strategy$name]] <- plan_df
    }

    if (save_objects) {
      saveRDS(qc_obj, file = file.path(strategy_dir, paste0(strategy$name, "_filtered.rds")))
    }
    filtered_objects[[strategy$name]] <- qc_obj
  }

  summary_df <- do.call(rbind, summary_rows)
  sample_retention_df <- do.call(rbind, sample_retention_rows)
  group_retention_df <- if (length(group_retention_rows)) do.call(rbind, group_retention_rows) else NULL
  cutoff_plan_df <- bind_rows_fill(cutoff_plan_rows)
  score_tables <- lapply(names(scoring_profile_list), function(profile_name) {
    score_benchmark_strategies(
      summary_df = summary_df,
      score_weights = scoring_profile_list[[profile_name]],
      scoring_profile = profile_name
    )
  })
  names(score_tables) <- names(scoring_profile_list)
  score_df <- do.call(rbind, score_tables)
  primary_score_df <- score_tables[[primary_profile]]
  recommendation_tables <- lapply(names(score_tables), function(profile_name) {
    rec <- recommend_benchmark_strategy(score_tables[[profile_name]])
    rec$scoring_profile <- profile_name
    rec[, c("scoring_profile", setdiff(colnames(rec), "scoring_profile")), drop = FALSE]
  })
  recommendations_df <- do.call(rbind, recommendation_tables)
  recommendation <- recommendations_df[recommendations_df$scoring_profile == primary_profile, , drop = FALSE]
  weights_used <- build_weights_table(scoring_profile_list)

  utils::write.csv(summary_df, file = file.path(output_dir, "benchmark_summary.csv"), row.names = FALSE)
  utils::write.csv(primary_score_df, file = file.path(output_dir, "benchmark_scores.csv"), row.names = FALSE)
  utils::write.csv(score_df, file = file.path(output_dir, "benchmark_normalized_scores.csv"), row.names = FALSE)
  utils::write.csv(weights_used, file = file.path(output_dir, "benchmark_score_weights.csv"), row.names = FALSE)
  utils::write.csv(recommendations_df, file = file.path(output_dir, "benchmark_recommended_strategies.csv"), row.names = FALSE)
  utils::write.csv(sample_retention_df, file = file.path(output_dir, "benchmark_sample_retention.csv"), row.names = FALSE)
  if (!is.null(group_retention_df)) {
    utils::write.csv(group_retention_df, file = file.path(output_dir, "benchmark_group_retention.csv"), row.names = FALSE)
  }
  if (!is.null(cutoff_plan_df)) {
    utils::write.csv(cutoff_plan_df, file = file.path(output_dir, "benchmark_cutoff_plan.csv"), row.names = FALSE)
  }

  write_benchmark_report(
    summary_df = summary_df,
    score_df = primary_score_df,
    recommendation = recommendation,
    score_weights = scoring_profile_list[[primary_profile]],
    weights_used = weights_used,
    recommendations_df = recommendations_df,
    sample_retention_df = sample_retention_df,
    group_retention_df = group_retention_df,
    output_dir = output_dir
  )

  if (plot) {
    benchmark_plots <- build_benchmark_plots(
      summary_df = summary_df,
      sample_retention_df = sample_retention_df,
      group_retention_df = group_retention_df,
      score_df = primary_score_df
    )
    for (plot_name in names(benchmark_plots)) {
      ggplot2::ggsave(
        filename = file.path(output_dir, paste0(plot_name, ".pdf")),
        plot = benchmark_plots[[plot_name]],
        width = if (identical(plot_name, "benchmark_sample_retention")) 10 else 8,
        height = if (identical(plot_name, "benchmark_sample_retention")) 5 else 4.8
      )
    }
  }

  invisible(list(
    summary = summary_df,
    scores = primary_score_df,
    recommendation = recommendation,
    raw_benchmark_table = summary_df,
    normalized_score_table = score_df,
    weights_used = weights_used,
    recommended_strategies = recommendations_df,
    sample_retention = sample_retention_df,
    group_retention = group_retention_df,
    cutoff_plan = cutoff_plan_df,
    strategy_objects = filtered_objects
  ))
}

normalize_benchmark_strategies <- function(strategies, group_by = NULL) {
  if (is.null(strategies)) {
    strategies <- default_benchmark_strategies(group_by = group_by)
  }

  lapply(strategies, function(strategy) {
    required <- c("name", "mode", "cutoff_strategy", "max_mito")
    missing_fields <- setdiff(required, names(strategy))
    if (length(missing_fields)) {
      stop(
        "Each strategy must define: ",
        paste(required, collapse = ", "),
        ". Missing: ",
        paste(missing_fields, collapse = ", "),
        call. = FALSE
      )
    }
    strategy$mode <- match.arg(strategy$mode, c("all", "split"))
    strategy$cutoff_strategy <- match.arg(strategy$cutoff_strategy, c("consensus", "strictest", "groupwise"))
    if (is.null(strategy$scdet_options)) {
      strategy$scdet_options <- list()
    }
    if (!is.list(strategy$scdet_options)) {
      stop("Each strategy 'scdet_options' entry must be a list.", call. = FALSE)
    }
    strategy
  })
}

bind_rows_fill <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(NULL)
  }
  all_columns <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing_columns <- setdiff(all_columns, colnames(row))
    for (column in missing_columns) {
      row[[column]] <- NA
    }
    row[, all_columns, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

default_benchmark_strategies <- function(group_by = NULL) {
  base <- list(
    list(name = "fixed_5", mode = "all", cutoff_strategy = "consensus", max_mito = 0.05),
    list(name = "fixed_10", mode = "all", cutoff_strategy = "consensus", max_mito = 0.10),
    list(name = "fixed_15", mode = "all", cutoff_strategy = "consensus", max_mito = 0.15),
    list(name = "fixed_20", mode = "all", cutoff_strategy = "consensus", max_mito = 0.20),
    list(
      name = "SCdetMito_largest_drop",
      mode = "all",
      cutoff_strategy = "consensus",
      max_mito = "SCdetMito",
      scdet_options = list(sample_cutoff_method = "largest_drop")
    ),
    list(
      name = "SCdetMito_first_significant_high",
      mode = "all",
      cutoff_strategy = "consensus",
      max_mito = "SCdetMito",
      scdet_options = list(sample_cutoff_method = "first_significant_high")
    ),
    list(
      name = "SCdetMito_reference_guided",
      mode = "all",
      cutoff_strategy = "consensus",
      max_mito = "SCdetMito",
      scdet_options = list(sample_cutoff_method = "reference_guided")
    ),
    list(name = "SCdetMito_strictest", mode = "all", cutoff_strategy = "strictest", max_mito = "SCdetMito")
  )

  if (!is.null(group_by)) {
    base[[length(base) + 1]] <- list(
      name = "SCdetMito_groupwise",
      mode = "all",
      cutoff_strategy = "groupwise",
      max_mito = "SCdetMito"
    )
  }

  base
}

summarize_benchmark_strategy <- function(original_obj,
                                         filtered_obj,
                                         strategy,
                                         sample_by,
                                         group_by,
                                         nFeature_RNA,
                                         nCount_RNA,
                                         mitoRatio,
                                         run_downstream,
                                         cluster_resolution,
                                         npcs,
                                         seed) {
  before_summary <- build_qc_summary(
    original_obj,
    nFeature_RNA = nFeature_RNA,
    nCount_RNA = nCount_RNA,
    mitoRatio = mitoRatio
  )
  after_summary <- build_qc_summary(
    filtered_obj,
    nFeature_RNA = nFeature_RNA,
    nCount_RNA = nCount_RNA,
    mitoRatio = mitoRatio
  )

  sample_retention <- summarize_retention_distribution(original_obj, filtered_obj, sample_by)
  group_retention <- if (is.null(group_by)) {
    empty_retention_distribution("group")
  } else {
    summarize_retention_distribution(original_obj, filtered_obj, group_by)
  }
  cutoff_stats <- extract_cutoff_stats(filtered_obj, strategy)
  detection_settings <- extract_detection_settings(filtered_obj)
  downstream_metrics <- if (isTRUE(run_downstream)) {
    compute_downstream_metrics(
      filtered_obj = filtered_obj,
      sample_by = sample_by,
      group_by = group_by,
      cluster_resolution = cluster_resolution,
      npcs = npcs,
      seed = seed
    )
  } else {
    empty_downstream_metrics("disabled")
  }

  data.frame(
    strategy = strategy$name,
    mode = strategy$mode,
    cutoff_strategy = strategy$cutoff_strategy,
    max_mito_input = as.character(strategy$max_mito),
    total_cells_before = before_summary$CellCounts,
    total_cells_after = after_summary$CellCounts,
    retention_rate = round(after_summary$CellCounts / before_summary$CellCounts, 4),
    retained_cell_fraction = round(after_summary$CellCounts / before_summary$CellCounts, 4),
    median_genes_before = before_summary$GenesCounts_median,
    median_genes_after = after_summary$GenesCounts_median,
    median_nFeature_RNA_after_qc = after_summary$GenesCounts_median,
    median_counts_before = before_summary$TranscriptCounts_median,
    median_counts_after = after_summary$TranscriptCounts_median,
    median_nCount_RNA_after_qc = after_summary$TranscriptCounts_median,
    median_mito_before = before_summary$MitoRatio_median,
    median_mito_after = after_summary$MitoRatio_median,
    median_mito_ratio_after_qc = after_summary$MitoRatio_median,
    sample_retention_min = sample_retention$retention_min,
    sample_retention_median = sample_retention$retention_median,
    sample_retention_max = sample_retention$retention_max,
    sample_retention_cv = sample_retention$retention_cv,
    sample_balance = retention_cv_to_balance(sample_retention$retention_cv),
    group_retention_min = group_retention$retention_min,
    group_retention_median = group_retention$retention_median,
    group_retention_max = group_retention$retention_max,
    group_retention_cv = group_retention$retention_cv,
    group_balance = retention_cv_to_balance(group_retention$retention_cv),
    applied_cutoff_min = cutoff_stats$AppliedCutoffMin,
    applied_cutoff_median = cutoff_stats$AppliedCutoffMedian,
    applied_cutoff_max = cutoff_stats$AppliedCutoffMax,
    loss_test = detection_settings$loss_test,
    p_adjust_method = detection_settings$p_adjust_method,
    alpha = detection_settings$alpha,
    sample_cutoff_method = detection_settings$sample_cutoff_method,
    cluster_count = downstream_metrics$cluster_count,
    sample_entropy_weighted = downstream_metrics$sample_entropy_weighted,
    group_entropy_weighted = downstream_metrics$group_entropy_weighted,
    dominant_sample_fraction = downstream_metrics$dominant_sample_fraction,
    dominant_group_fraction = downstream_metrics$dominant_group_fraction,
    downstream_status = downstream_metrics$status,
    stringsAsFactors = FALSE
  )
}

build_retention_table <- function(original_obj,
                                  filtered_obj,
                                  column,
                                  level_type,
                                  strategy_name) {
  before_counts <- table(as.character(original_obj@meta.data[[column]]))
  after_counts <- table(as.character(filtered_obj@meta.data[[column]]))
  levels <- union(names(before_counts), names(after_counts))
  before_vec <- as.numeric(before_counts[levels])
  after_vec <- as.numeric(after_counts[levels])
  before_vec[is.na(before_vec)] <- 0
  after_vec[is.na(after_vec)] <- 0

  data.frame(
    strategy = strategy_name,
    level_type = level_type,
    level_id = levels,
    before_cells = before_vec,
    after_cells = after_vec,
    retention_rate = round(after_vec / before_vec, 4),
    stringsAsFactors = FALSE
  )
}

summarize_retention_distribution <- function(original_obj, filtered_obj, column) {
  retention_table <- build_retention_table(
    original_obj = original_obj,
    filtered_obj = filtered_obj,
    column = column,
    level_type = column,
    strategy_name = "temp"
  )
  rates <- retention_table$retention_rate[is.finite(retention_table$retention_rate)]
  if (!length(rates)) {
    return(empty_retention_distribution(column))
  }

  data.frame(
    retention_min = min(rates, na.rm = TRUE),
    retention_median = stats::median(rates, na.rm = TRUE),
    retention_max = max(rates, na.rm = TRUE),
    retention_cv = if (length(rates) <= 1) {
      0
    } else if (mean(rates, na.rm = TRUE) == 0) {
      NA_real_
    } else {
      stats::sd(rates, na.rm = TRUE) / mean(rates, na.rm = TRUE)
    },
    stringsAsFactors = FALSE
  )
}

empty_retention_distribution <- function(prefix) {
  data.frame(
    retention_min = NA_real_,
    retention_median = NA_real_,
    retention_max = NA_real_,
    retention_cv = NA_real_,
    stringsAsFactors = FALSE
  )
}

retention_cv_to_balance <- function(retention_cv) {
  ifelse(is.finite(retention_cv), round(1 / (1 + retention_cv), 4), NA_real_)
}

extract_cutoff_stats <- function(filtered_obj, strategy) {
  cutoff_plan <- filtered_obj@misc$SCdetMito_QC
  if (is.null(cutoff_plan$sample_plan)) {
    fixed_cutoff <- suppressWarnings(as.numeric(strategy$max_mito))
    return(data.frame(
      AppliedCutoffMin = fixed_cutoff,
      AppliedCutoffMedian = fixed_cutoff,
      AppliedCutoffMax = fixed_cutoff,
      stringsAsFactors = FALSE
    ))
  }
  summarize_cutoff_vector(cutoff_plan$sample_plan$applied_cutoff)
}

extract_detection_settings <- function(filtered_obj) {
  detection_settings <- filtered_obj@misc$SCdetMito_QC$detection$settings
  if (is.null(detection_settings)) {
    return(data.frame(
      loss_test = NA_character_,
      p_adjust_method = NA_character_,
      alpha = NA_real_,
      sample_cutoff_method = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    loss_test = detection_settings$loss_test %||% NA_character_,
    p_adjust_method = detection_settings$p_adjust_method %||% NA_character_,
    alpha = detection_settings$alpha %||% NA_real_,
    sample_cutoff_method = detection_settings$sample_cutoff_method %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

compute_downstream_metrics <- function(filtered_obj,
                                       sample_by,
                                       group_by,
                                       cluster_resolution,
                                       npcs,
                                       seed) {
  if (ncol(filtered_obj) < 20) {
    return(empty_downstream_metrics("too_few_cells"))
  }

  result <- tryCatch(
    {
      set.seed(seed)
      seu <- Seurat::NormalizeData(filtered_obj, verbose = FALSE)
      seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      seu <- Seurat::ScaleData(seu, verbose = FALSE)
      max_pc <- min(as.integer(npcs), ncol(seu) - 1L)
      if (max_pc < 5L) {
        return(empty_downstream_metrics("too_few_pcs"))
      }
      seu <- Seurat::RunPCA(seu, npcs = max_pc, verbose = FALSE)
      dims <- seq_len(max_pc)
      seu <- Seurat::FindNeighbors(seu, dims = dims, verbose = FALSE)
      seu <- Seurat::FindClusters(seu, resolution = cluster_resolution, verbose = FALSE)

      cluster_labels <- as.character(seu$seurat_clusters)
      data.frame(
        cluster_count = length(unique(cluster_labels)),
        sample_entropy_weighted = compute_weighted_cluster_entropy(cluster_labels, seu@meta.data[[sample_by]]),
        group_entropy_weighted = if (is.null(group_by)) {
          NA_real_
        } else {
          compute_weighted_cluster_entropy(cluster_labels, seu@meta.data[[group_by]])
        },
        dominant_sample_fraction = compute_weighted_dominant_fraction(cluster_labels, seu@meta.data[[sample_by]]),
        dominant_group_fraction = if (is.null(group_by)) {
          NA_real_
        } else {
          compute_weighted_dominant_fraction(cluster_labels, seu@meta.data[[group_by]])
        },
        status = "ok",
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      empty_downstream_metrics("failed")
    }
  )

  result
}

empty_downstream_metrics <- function(status) {
  data.frame(
    cluster_count = NA_real_,
    sample_entropy_weighted = NA_real_,
    group_entropy_weighted = NA_real_,
    dominant_sample_fraction = NA_real_,
    dominant_group_fraction = NA_real_,
    status = status,
    stringsAsFactors = FALSE
  )
}

compute_weighted_cluster_entropy <- function(clusters, labels) {
  cluster_ids <- unique(as.character(clusters))
  cluster_sizes <- table(clusters)
  entropy_values <- vapply(cluster_ids, function(cluster_id) {
    label_counts <- table(as.character(labels[clusters == cluster_id]))
    proportions <- as.numeric(label_counts) / sum(label_counts)
    if (length(proportions) <= 1) {
      return(0)
    }
    -sum(proportions * log(proportions)) / log(length(proportions))
  }, numeric(1))
  stats::weighted.mean(entropy_values, w = as.numeric(cluster_sizes[cluster_ids]), na.rm = TRUE)
}

compute_weighted_dominant_fraction <- function(clusters, labels) {
  cluster_ids <- unique(as.character(clusters))
  cluster_sizes <- table(clusters)
  dominant_values <- vapply(cluster_ids, function(cluster_id) {
    label_counts <- table(as.character(labels[clusters == cluster_id]))
    max(label_counts) / sum(label_counts)
  }, numeric(1))
  stats::weighted.mean(dominant_values, w = as.numeric(cluster_sizes[cluster_ids]), na.rm = TRUE)
}

build_benchmark_plots <- function(summary_df,
                                  sample_retention_df,
                                  group_retention_df = NULL,
                                  score_df = NULL) {
  plots <- list()
  cells_df <- summary_df[order(summary_df$total_cells_after), , drop = FALSE]
  cells_df$strategy_order <- factor(cells_df$strategy, levels = cells_df$strategy)
  mito_df <- summary_df[order(summary_df$median_mito_after), , drop = FALSE]
  mito_df$strategy_order <- factor(mito_df$strategy, levels = mito_df$strategy)
  palette_values <- scdetmito_palette(length(unique(summary_df$cutoff_strategy)))

  plots$benchmark_cells_retained <- ggplot2::ggplot(
    cells_df,
    ggplot2::aes(x = strategy_order, y = total_cells_after, fill = cutoff_strategy)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = palette_values) +
    theme_scdetmito() +
    ggplot2::labs(x = "Strategy", y = "Retained cells", fill = "Cutoff strategy")

  plots$benchmark_median_mito <- ggplot2::ggplot(
    mito_df,
    ggplot2::aes(x = strategy_order, y = median_mito_after, fill = cutoff_strategy)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = palette_values) +
    theme_scdetmito() +
    ggplot2::labs(x = "Strategy", y = "Median mitoRatio after QC", fill = "Cutoff strategy")

  plots$benchmark_sample_retention <- ggplot2::ggplot(
    sample_retention_df,
    ggplot2::aes(x = level_id, y = strategy, fill = retention_rate)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "#f6f7fb", high = "#3b6ea8") +
    theme_scdetmito() +
    ggplot2::labs(x = "Sample", y = "Strategy", fill = "Retention")

  if (!is.null(group_retention_df)) {
    plots$benchmark_group_retention <- ggplot2::ggplot(
      group_retention_df,
      ggplot2::aes(x = level_id, y = strategy, fill = retention_rate)
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low = "#f6f7fb", high = "#2f7d4a") +
      theme_scdetmito() +
      ggplot2::labs(x = "Group", y = "Strategy", fill = "Retention")
  }

  if (!is.null(score_df) && nrow(score_df)) {
    score_plot_df <- score_df[order(score_df$overall_score), , drop = FALSE]
    score_plot_df$strategy_order <- factor(score_plot_df$strategy, levels = score_plot_df$strategy)
    plots$SCQCbenchmark_strategy_comparison <- ggplot2::ggplot(
      score_plot_df,
      ggplot2::aes(x = strategy_order, y = overall_score, fill = cutoff_strategy)
    ) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = palette_values) +
      theme_scdetmito() +
      ggplot2::labs(
        x = "Strategy",
        y = "Composite score under selected profile",
        fill = "Cutoff strategy"
      )
  }

  plots
}

write_benchmark_report <- function(summary_df,
                                   score_df,
                                   recommendation,
                                   score_weights,
                                   weights_used = NULL,
                                   recommendations_df = NULL,
                                   sample_retention_df,
                                   group_retention_df = NULL,
                                   output_dir = ".") {
  highest_retention <- summary_df$strategy[which.max(summary_df$retention_rate)]
  lowest_sample_cv <- summary_df$strategy[which.min(summary_df$sample_retention_cv)]
  top_scoring_strategy <- score_df$strategy[which.max(score_df$overall_score)]
  recommended_adaptive_strategy <- recommendation$recommended_adaptive_strategy %||%
    recommendation$adaptive_strategy

  report_lines <- c(
    "# QC Benchmark Report",
    "",
    paste0("- Highest cell retention: `", highest_retention, "` (", max(summary_df$retention_rate), ")"),
    paste0("- Lowest sample-retention CV: `", lowest_sample_cv, "` (", round(min(summary_df$sample_retention_cv, na.rm = TRUE), 4), ")"),
    paste0("- Top-scoring strategy under selected profile: `", top_scoring_strategy, "`"),
    if (!all(is.na(summary_df$group_retention_cv))) {
      paste0(
        "- Lowest group-retention CV: `",
        summary_df$strategy[which.min(summary_df$group_retention_cv)],
        "` (",
        round(min(summary_df$group_retention_cv, na.rm = TRUE), 4),
        ")"
      )
    },
    if (!is.na(recommended_adaptive_strategy)) {
      paste0("- Recommended adaptive strategy under selected profile: `", recommended_adaptive_strategy, "`")
    },
    paste0(
      "- Score weights: retention=",
      round(score_weights[["retention"]], 3),
      ", mito=",
      round(score_weights[["mito"]], 3),
      ", sample_balance=",
      round(score_weights[["sample_balance"]], 3),
      ", group_balance=",
      round(score_weights[["group_balance"]], 3)
    ),
    "",
    "## Recommendation",
    "",
    if (!is.na(recommended_adaptive_strategy)) {
      paste0(
        "`",
        recommended_adaptive_strategy,
        "` has the highest adaptive composite score under the selected scoring profile for this dataset."
      )
    },
    if (!is.na(recommended_adaptive_strategy)) {
      paste0(
        "Key metrics: retention = ",
        recommendation$retention_rate,
        ", median mito after QC = ",
        recommendation$median_mito_after,
        ", sample CV = ",
        recommendation$sample_retention_cv,
        if (!is.na(recommendation$group_retention_cv)) {
          paste0(", group CV = ", recommendation$group_retention_cv)
        } else {
          ""
        },
        "."
      )
    },
    "",
    "## Strategy Summary",
    ""
  )

  summary_table <- paste(utils::capture.output(print(summary_df, row.names = FALSE)), collapse = "\n")
  report_lines <- c(report_lines, "```text", summary_table, "```")

  report_lines <- c(report_lines, "", "## Strategy Scores", "")
  score_table <- paste(utils::capture.output(print(score_df, row.names = FALSE)), collapse = "\n")
  report_lines <- c(report_lines, "```text", score_table, "```")

  if (!is.null(recommendations_df)) {
    report_lines <- c(report_lines, "", "## Scoring Profiles", "")
    recommendation_table <- paste(utils::capture.output(print(recommendations_df, row.names = FALSE)), collapse = "\n")
    report_lines <- c(report_lines, "```text", recommendation_table, "```")
  }

  if (!is.null(weights_used)) {
    report_lines <- c(report_lines, "", "## Score Weights", "")
    weights_table <- paste(utils::capture.output(print(weights_used, row.names = FALSE)), collapse = "\n")
    report_lines <- c(report_lines, "```text", weights_table, "```")
  }

  if (!is.null(group_retention_df)) {
    report_lines <- c(report_lines, "", "## Group Retention", "")
    group_table <- paste(utils::capture.output(print(group_retention_df, row.names = FALSE)), collapse = "\n")
    report_lines <- c(report_lines, "```text", group_table, "```")
  }

  writeLines(report_lines, con = file.path(output_dir, "benchmark_report.md"))
}

score_benchmark_strategies <- function(summary_df, score_weights, scoring_profile = "balanced") {
  score_df <- summary_df
  score_df$scoring_profile <- scoring_profile
  score_df$retention_score <- scale_benchmark_metric(score_df$retention_rate, higher_better = TRUE)
  score_df$mito_score <- scale_benchmark_metric(score_df$median_mito_after, higher_better = FALSE)
  score_df$sample_balance_score <- scale_benchmark_metric(score_df$sample_retention_cv, higher_better = FALSE)
  score_df$group_balance_score <- scale_benchmark_metric(
    ifelse(is.na(score_df$group_retention_cv), score_df$sample_retention_cv, score_df$group_retention_cv),
    higher_better = FALSE
  )
  score_matrix <- cbind(
    retention = score_df$retention_score,
    mito = score_df$mito_score,
    sample_balance = score_df$sample_balance_score,
    group_balance = score_df$group_balance_score
  )
  score_df$overall_score <- round(apply(score_matrix, 1, function(row_scores) {
    valid <- is.finite(row_scores)
    if (!any(valid)) {
      return(NA_real_)
    }
    sum(row_scores[valid] * score_weights[names(row_scores)[valid]]) /
      sum(score_weights[names(row_scores)[valid]])
  }), 4)
  score_df[order(-score_df$overall_score, -score_df$retention_rate), , drop = FALSE]
}

scale_benchmark_metric <- function(values, higher_better = TRUE) {
  clean_values <- values
  clean_values[!is.finite(clean_values)] <- NA_real_

  if (all(is.na(clean_values))) {
    return(rep(NA_real_, length(values)))
  }

  range_values <- range(clean_values, na.rm = TRUE)
  if (diff(range_values) == 0) {
    scaled <- rep(1, length(values))
  } else if (isTRUE(higher_better)) {
    scaled <- (clean_values - range_values[1]) / diff(range_values)
  } else {
    scaled <- (range_values[2] - clean_values) / diff(range_values)
  }

  scaled[is.na(clean_values)] <- NA_real_
  round(scaled, 4)
}

recommend_benchmark_strategy <- function(score_df) {
  adaptive_df <- score_df[grepl("^(scdet_|SCdetMito_)", as.character(score_df$strategy)), , drop = FALSE]
  if (!nrow(adaptive_df)) {
    return(data.frame(
      recommended_adaptive_strategy = NA_character_,
      recommended_strategy = score_df$strategy[1],
      top_scoring_strategy = score_df$strategy[1],
      adaptive_strategy = NA_character_,
      overall_strategy = score_df$strategy[1],
      retention_rate = NA_real_,
      median_mito_after = NA_real_,
      sample_retention_cv = NA_real_,
      group_retention_cv = NA_real_,
      overall_score = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  adaptive_scores <- adaptive_df$overall_score
  valid_idx <- which(is.finite(adaptive_scores))
  recommended_adaptive <- if (length(valid_idx)) {
    adaptive_df[valid_idx[which.max(adaptive_scores[valid_idx])], , drop = FALSE]
  } else {
    adaptive_df[1, , drop = FALSE]
  }
  data.frame(
    recommended_adaptive_strategy = recommended_adaptive$strategy,
    recommended_strategy = score_df$strategy[1],
    top_scoring_strategy = score_df$strategy[1],
    adaptive_strategy = recommended_adaptive$strategy,
    overall_strategy = score_df$strategy[1],
    retention_rate = recommended_adaptive$retention_rate,
    median_mito_after = recommended_adaptive$median_mito_after,
    sample_retention_cv = recommended_adaptive$sample_retention_cv,
    group_retention_cv = recommended_adaptive$group_retention_cv,
    overall_score = recommended_adaptive$overall_score,
    stringsAsFactors = FALSE
  )
}

normalize_score_weights <- function(score_weights) {
  required_names <- c("retention", "mito", "sample_balance", "group_balance")
  if (is.list(score_weights)) {
    score_weights <- unlist(score_weights, use.names = TRUE)
  }
  if (!is.numeric(score_weights) || is.null(names(score_weights))) {
    stop(
      "'score_weights' must be a named numeric vector or list with: ",
      paste(required_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!all(required_names %in% names(score_weights))) {
    stop(
      "'score_weights' is missing required entries: ",
      paste(setdiff(required_names, names(score_weights)), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  normalized <- as.numeric(score_weights[required_names])
  if (any(!is.finite(normalized)) || any(normalized < 0)) {
    stop("'score_weights' values must be finite and non-negative.", call. = FALSE)
  }
  total_weight <- sum(normalized)
  if (total_weight <= 0) {
    stop("'score_weights' must sum to a positive value.", call. = FALSE)
  }
  normalized <- normalized / total_weight
  names(normalized) <- required_names
  normalized
}

resolve_scoring_profiles <- function(scoring_profiles, score_weights = NULL) {
  built_in <- list(
    balanced = c(retention = 0.35, mito = 0.25, sample_balance = 0.20, group_balance = 0.20),
    retention_focused = c(retention = 0.60, mito = 0.15, sample_balance = 0.15, group_balance = 0.10),
    stringent_mito = c(retention = 0.20, mito = 0.55, sample_balance = 0.15, group_balance = 0.10),
    balance_focused = c(retention = 0.20, mito = 0.15, sample_balance = 0.35, group_balance = 0.30)
  )

  if (is.null(scoring_profiles)) {
    scoring_profiles <- names(built_in)
  }
  scoring_profiles <- unique(match.arg(
    scoring_profiles,
    choices = names(built_in),
    several.ok = TRUE
  ))
  profiles <- lapply(built_in[scoring_profiles], normalize_score_weights)
  if (!is.null(score_weights)) {
    profiles <- c(custom = list(normalize_score_weights(score_weights)), profiles)
  }
  profiles
}

build_weights_table <- function(scoring_profile_list) {
  rows <- lapply(names(scoring_profile_list), function(profile_name) {
    weights <- scoring_profile_list[[profile_name]]
    data.frame(
      scoring_profile = profile_name,
      retention = weights[["retention"]],
      mito = weights[["mito"]],
      sample_balance = weights[["sample_balance"]],
      group_balance = weights[["group_balance"]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
