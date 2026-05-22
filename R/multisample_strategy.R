# Internal utilities for group-aware multi-sample QC planning.

normalize_group_by <- function(group_by) {
  if (is.null(group_by) || identical(group_by, "") || is.na(group_by)) {
    return(NULL)
  }
  as.character(group_by)
}

summarize_cutoff_vector <- function(applied_cutoffs) {
  data.frame(
    AppliedCutoffMin = min(applied_cutoffs, na.rm = TRUE),
    AppliedCutoffMedian = stats::median(applied_cutoffs, na.rm = TRUE),
    AppliedCutoffMax = max(applied_cutoffs, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

add_cutoff_metrics <- function(summary_df, applied_cutoffs, cutoff_strategy) {
  cutoff_summary <- summarize_cutoff_vector(applied_cutoffs)
  summary_df$CutoffStrategy <- cutoff_strategy
  cbind(summary_df, cutoff_summary)
}

build_sample_group_map <- function(seurat_obj, sample_by, group_by = NULL) {
  metadata <- seurat_obj@meta.data
  sample_ids <- unique(as.character(metadata[[sample_by]]))
  sample_map <- data.frame(
    sample_id = sample_ids,
    raw_cells = vapply(sample_ids, function(sample_id) {
      sum(metadata[[sample_by]] == sample_id, na.rm = TRUE)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )

  if (is.null(group_by)) {
    sample_map$group_id <- NA_character_
    return(sample_map)
  }

  sample_groups <- lapply(sample_ids, function(sample_id) {
    unique(as.character(stats::na.omit(metadata[metadata[[sample_by]] == sample_id, group_by])))
  })
  invalid_mapping <- lengths(sample_groups) != 1
  if (any(invalid_mapping)) {
    stop(
      "Each sample in '", sample_by, "' must map to exactly one value in '", group_by, "'.",
      call. = FALSE
    )
  }

  sample_map$group_id <- vapply(sample_groups, `[`, character(1), 1)
  sample_map
}

aggregate_cutoff_value <- function(cutoffs, cutoff_quantile = 0.25) {
  as.numeric(stats::quantile(
    cutoffs,
    probs = cutoff_quantile,
    na.rm = TRUE,
    names = FALSE,
    type = 8
  ))
}

resolve_supported_cutoff <- function(cutoffs,
                                     min_support_fraction = 0.5,
                                     n_total_samples = length(cutoffs),
                                     strategy = "sample_supported_global_cutoff",
                                     group = NA_character_) {
  raw_cutoffs <- as.numeric(cutoffs)
  cutoffs <- raw_cutoffs[is.finite(raw_cutoffs)]
  if (!is.numeric(min_support_fraction) ||
    length(min_support_fraction) != 1 ||
    is.na(min_support_fraction) ||
    min_support_fraction <= 0 ||
    min_support_fraction > 1) {
      stop("'min_support_fraction' must be a numeric scalar in the interval (0, 1].", call. = FALSE)
  }

  required_support <- max(1L, ceiling(n_total_samples * min_support_fraction))
  if (!length(cutoffs)) {
    return(data.frame(
      global_cutoff = NA_real_,
      cutoff = NA_real_,
      support_fraction = 0,
      n_supporting_samples = 0L,
      support_count = 0L,
      n_total_samples = n_total_samples,
      sample_count = n_total_samples,
      required_support = required_support,
      strategy = strategy,
      group = group,
      stringsAsFactors = FALSE
    ))
  }
  candidate_cutoffs <- sort(unique(cutoffs))
  support_counts <- vapply(candidate_cutoffs, function(candidate) {
    sum(cutoffs >= candidate, na.rm = TRUE)
  }, integer(1))
  supported <- support_counts >= required_support
  selected <- if (any(supported)) {
    max(candidate_cutoffs[supported], na.rm = TRUE)
  } else {
    min(candidate_cutoffs, na.rm = TRUE)
  }
  selected_support <- support_counts[match(selected, candidate_cutoffs)]

  data.frame(
    global_cutoff = selected,
    cutoff = selected,
    support_fraction = selected_support / n_total_samples,
    n_supporting_samples = selected_support,
    support_count = selected_support,
    n_total_samples = n_total_samples,
    sample_count = n_total_samples,
    required_support = required_support,
    strategy = strategy,
    group = group,
    stringsAsFactors = FALSE
  )
}

summarize_cutoff_support <- function(cutoffs,
                                     cutoff,
                                     min_support_fraction = 0.5,
                                     n_total_samples = length(cutoffs),
                                     strategy = "sample_supported_global_cutoff",
                                     group = NA_character_) {
  raw_cutoffs <- as.numeric(cutoffs)
  cutoffs <- raw_cutoffs[is.finite(raw_cutoffs)]
  if (!is.numeric(cutoff) ||
    length(cutoff) != 1 ||
    !is.finite(cutoff)) {
    return(data.frame(
      global_cutoff = NA_real_,
      cutoff = NA_real_,
      support_fraction = 0,
      n_supporting_samples = 0L,
      support_count = 0L,
      n_total_samples = n_total_samples,
      sample_count = n_total_samples,
      required_support = max(1L, ceiling(n_total_samples * min_support_fraction)),
      strategy = strategy,
      group = group,
      stringsAsFactors = FALSE
    ))
  }
  if (!is.numeric(min_support_fraction) ||
    length(min_support_fraction) != 1 ||
    is.na(min_support_fraction) ||
    min_support_fraction <= 0 ||
    min_support_fraction > 1) {
      stop("'min_support_fraction' must be a numeric scalar in the interval (0, 1].", call. = FALSE)
  }

  required_support <- max(1L, ceiling(n_total_samples * min_support_fraction))
  support_count <- sum(cutoffs >= cutoff, na.rm = TRUE)
  data.frame(
    global_cutoff = cutoff,
    cutoff = cutoff,
    support_fraction = support_count / n_total_samples,
    n_supporting_samples = support_count,
    support_count = support_count,
    n_total_samples = n_total_samples,
    sample_count = n_total_samples,
    required_support = required_support,
    strategy = strategy,
    group = group,
    stringsAsFactors = FALSE
  )
}

build_multisample_cutoff_plan <- function(seurat_obj,
                                          sample_by,
                                          group_by = NULL,
                                          mitoRatio = "mitoRatio",
                                          max_mito = "SCdetMito",
                                          cutoff_strategy = c("consensus", "strictest", "groupwise"),
                                          cutoff_quantile = 0.25,
                                          cutoff_support_fraction = 0.5,
                                          scdet_options = list(),
                                          use_recommended_cutoff = TRUE,
                                          table_out = TRUE,
                                          plot = TRUE,
                                          output_dir = ".") {
  cutoff_strategy <- match.arg(cutoff_strategy)
  group_by <- normalize_group_by(group_by)

  if (!is.numeric(cutoff_quantile) ||
    length(cutoff_quantile) != 1 ||
    is.na(cutoff_quantile) ||
    cutoff_quantile <= 0 ||
    cutoff_quantile > 1) {
    stop("'cutoff_quantile' must be a numeric scalar in the interval (0, 1].", call. = FALSE)
  }
  if (!is.list(scdet_options)) {
    stop("'scdet_options' must be a list.", call. = FALSE)
  }
  if (!is.numeric(cutoff_support_fraction) ||
    length(cutoff_support_fraction) != 1 ||
    is.na(cutoff_support_fraction) ||
    cutoff_support_fraction <= 0 ||
    cutoff_support_fraction > 1) {
    stop("'cutoff_support_fraction' must be a numeric scalar in the interval (0, 1].", call. = FALSE)
  }

  sample_plan <- build_sample_group_map(
    seurat_obj = seurat_obj,
    sample_by = sample_by,
    group_by = group_by
  )

  detection <- NULL
  if (identical(max_mito, "SCdetMito")) {
    detection_args <- utils::modifyList(list(
      seurat_obj = seurat_obj,
      mitoRatio = mitoRatio,
      by = sample_by,
      table_out = table_out,
      plot = plot,
      output_dir = output_dir,
      return_details = TRUE
    ), scdet_options)
    detection <- do.call(
      SCdetMito,
      detection_args
    )
    detected_summary <- detection$sample_cutoff_summary
    colnames(detected_summary)[colnames(detected_summary) %in% c(sample_by, "sample")] <- "sample_id"
    sample_plan <- merge(
      sample_plan,
      detected_summary,
      by = "sample_id",
      all.x = TRUE,
      sort = FALSE
    )
    sample_plan <- sample_plan[match(unique(as.character(seurat_obj@meta.data[[sample_by]])), sample_plan$sample_id), ]
    if (anyNA(sample_plan$detected_cutoff)) {
      stop("Failed to assign a detected cutoff to every sample.", call. = FALSE)
    }
    cutoff_column <- if (isTRUE(use_recommended_cutoff) &&
      "recommended_cutoff" %in% colnames(sample_plan) &&
      any(is.finite(sample_plan$recommended_cutoff))) {
        "recommended_cutoff"
      } else {
        "selected_cutoff"
      }
    if (!cutoff_column %in% colnames(sample_plan)) {
      cutoff_column <- "detected_cutoff"
    }
    sample_plan$sample_cutoff <- sample_plan[[cutoff_column]]
    sample_plan$cutoff_applied_source <- if (identical(cutoff_column, "recommended_cutoff")) {
      "recommended_cutoff"
    } else {
      "selected_cutoff"
    }
  } else {
    numeric_cutoff <- normalize_mito_cutoff_value(
      max_mito,
      name = "max_mito",
      allow_scdet = FALSE
    )
    sample_plan$detected_cutoff <- numeric_cutoff
    sample_plan$selected_cutoff <- numeric_cutoff
    sample_plan$recommended_cutoff <- numeric_cutoff
    sample_plan$recommended_method <- "user_defined"
    sample_plan$recommendation_level <- "standard"
    sample_plan$recommendation_source <- "user_defined"
    sample_plan$sample_cutoff <- numeric_cutoff
    sample_plan$cutoff_applied_source <- "user_defined"
    sample_plan$cutoff_source <- "user_defined"
    sample_plan$significant_cutoff_count <- 0
    sample_plan$significant_interval_count <- 0
    sample_plan$fallback_used <- FALSE
    sample_plan$fallback_method <- NA_character_
    sample_plan$fallback_quantile <- NA_real_
    sample_plan$cutoff_confidence <- "user_defined"
  }

  if (!identical(cutoff_strategy, "groupwise")) {
    if (identical(cutoff_strategy, "strictest")) {
      final_cutoff <- min(sample_plan$sample_cutoff, na.rm = TRUE)
      supported_cutoff <- summarize_cutoff_support(
        sample_plan$sample_cutoff,
        final_cutoff,
        cutoff_support_fraction
      )
    } else {
      supported_cutoff <- resolve_supported_cutoff(
        sample_plan$sample_cutoff,
        cutoff_support_fraction
      )
      final_cutoff <- supported_cutoff$cutoff
    }
    sample_plan$group_cutoff <- final_cutoff
    sample_plan$applied_cutoff <- final_cutoff
    aggregation_level <- "global"
    supported_cutoff <- summarize_cutoff_support(
      sample_plan$sample_cutoff,
      final_cutoff,
      cutoff_support_fraction
    )
    group_plan <- data.frame(
      group_id = "all_samples",
      group_cutoff = if (is.na(final_cutoff)) NA_real_ else final_cutoff,
      sample_count = nrow(sample_plan),
      support_count = supported_cutoff$support_count,
      support_fraction = supported_cutoff$support_fraction,
      required_support = supported_cutoff$required_support,
      stringsAsFactors = FALSE
    )
  } else if (is.null(group_by)) {
    supported_cutoff <- resolve_supported_cutoff(sample_plan$sample_cutoff, cutoff_support_fraction)
    final_cutoff <- supported_cutoff$cutoff
    sample_plan$group_cutoff <- sample_plan$sample_cutoff
    sample_plan$applied_cutoff <- sample_plan$sample_cutoff
    aggregation_level <- "sample"
    group_plan <- data.frame(
      group_id = "all_samples",
      group_cutoff = if (is.na(final_cutoff)) NA_real_ else final_cutoff,
      sample_count = nrow(sample_plan),
      support_count = supported_cutoff$support_count,
      support_fraction = supported_cutoff$support_fraction,
      required_support = supported_cutoff$required_support,
      stringsAsFactors = FALSE
    )
  } else {
    group_cutoffs <- tapply(sample_plan$sample_cutoff, sample_plan$group_id, function(cutoffs) {
      resolve_supported_cutoff(cutoffs, cutoff_support_fraction)$cutoff
    })
    group_support <- do.call(
      rbind,
      lapply(names(group_cutoffs), function(group_id) {
        group_sample_cutoffs <- sample_plan$sample_cutoff[sample_plan$group_id == group_id]
        support <- summarize_cutoff_support(
          group_sample_cutoffs,
          as.numeric(group_cutoffs[[group_id]]),
          cutoff_support_fraction
        )
        data.frame(
          group_id = group_id,
          support_count = support$support_count,
          support_fraction = support$support_fraction,
          required_support = support$required_support,
          stringsAsFactors = FALSE
        )
      })
    )
    group_plan <- data.frame(
      group_id = names(group_cutoffs),
      group_cutoff = as.numeric(group_cutoffs),
      sample_count = vapply(names(group_cutoffs), function(group_id) {
        sum(sample_plan$group_id == group_id, na.rm = TRUE)
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
    group_plan <- merge(group_plan, group_support, by = "group_id", all.x = TRUE, sort = FALSE)
    group_plan <- group_plan[match(names(group_cutoffs), group_plan$group_id), ]

    sample_plan$group_cutoff <- group_plan$group_cutoff[match(sample_plan$group_id, group_plan$group_id)]

    final_cutoff <- max(group_plan$group_cutoff, na.rm = TRUE)
    sample_plan$applied_cutoff <- sample_plan$group_cutoff
    aggregation_level <- "group"
  }

  sample_plan$cutoff_strategy <- cutoff_strategy
  sample_plan$aggregation_level <- aggregation_level
  sample_plan$sample_by <- sample_by
  sample_plan$group_by <- if (is.null(group_by)) NA_character_ else group_by
  sample_plan$supports_final_cutoff <- sample_plan$sample_cutoff >= final_cutoff

  group_plan$global_cutoff <- group_plan$group_cutoff
  group_plan$n_supporting_samples <- group_plan$support_count
  group_plan$n_total_samples <- group_plan$sample_count
  group_plan$strategy <- cutoff_strategy
  group_plan$group <- group_plan$group_id

  if (anyNA(sample_plan$applied_cutoff)) {
    stop("Failed to derive applied mitochondrial cutoffs for all samples.", call. = FALSE)
  }

  list(
    sample_plan = sample_plan,
    group_plan = group_plan,
    final_cutoff = final_cutoff,
    global_cutoff = final_cutoff,
    sample_supported_global_cutoff = group_plan,
    cutoff_strategy = cutoff_strategy,
    strategy = cutoff_strategy,
    sample_by = sample_by,
    group_by = group_by,
    detection = detection
  )
}

build_cell_level_cutoff_vector <- function(seurat_obj, sample_by, sample_plan) {
  metadata <- seurat_obj@meta.data
  matched_cutoffs <- sample_plan$applied_cutoff[
    match(as.character(metadata[[sample_by]]), sample_plan$sample_id)
  ]
  if (anyNA(matched_cutoffs)) {
    stop("Failed to map applied cutoffs to all cells.", call. = FALSE)
  }
  names(matched_cutoffs) <- rownames(metadata)
  matched_cutoffs
}
