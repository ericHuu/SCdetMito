make_threshold_safety_seurat <- function(n_cells = 30L) {
  set.seed(20260909)
  counts <- matrix(
    rpois(40 * n_cells, lambda = 5),
    nrow = 40,
    ncol = n_cells
  )
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("SafetyCell", seq_len(ncol(counts)))
  seu <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = counts,
    min.cells = 0,
    min.features = 0
  ))
  seu$sample <- "A"
  seu
}

test_that("SCdetMito rejects invalid cutoff search ranges", {
  seu <- make_threshold_safety_seurat()
  seu$mitoRatio <- 0.05

  expect_error(
    SCdetMito(seu, sample_col = "sample", max_cut = 100, plot = FALSE),
    "max_cut"
  )
  expect_error(
    SCdetMito(seu, sample_col = "sample", min_cut = -0.01, plot = FALSE),
    "min_cut"
  )
  expect_error(
    SCdetMito(seu, sample_col = "sample", min_cut = c(0.01, 0.02), plot = FALSE),
    "min_cut"
  )
})

test_that("retained-cell profiles use the same inclusive cutoff as filtering", {
  seu <- make_threshold_safety_seurat(n_cells = 3L)
  seu$mitoRatio <- c(0.10, 0.10, 0.15)

  result <- suppressWarnings(SCdetMito(
    seu,
    sample_col = "sample",
    min_cut = 0.10,
    max_cut = 0.20,
    bin_width = 0.10,
    min_drop_cells = 100,
    min_cells_after = 1,
    fallback_method = "none",
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  ))
  at_point_one <- result$retained_cell_profile[
    result$retained_cell_profile$cutoff == 0.10,
    ,
    drop = FALSE
  ]

  expect_equal(at_point_one$retained_cells, 2L)
  expect_equal(result$settings$retained_cell_profile_formula, "R_s(c) = sum(m_i <= c)")
})

test_that("existing mitochondrial metadata supports explicit scales and rejects mixtures", {
  seu <- make_threshold_safety_seurat(n_cells = 20L)
  seu$mitoPercent <- seq(0.2, 0.8, length.out = 20)

  percent_result <- ensure_mito_ratio(
    seu,
    mito_col = "mitoPercent",
    mito_scale = "percent",
    verbose = FALSE
  )
  expect_equal(percent_result$mitoPercent, seu$mitoPercent / 100)
  expect_equal(percent_result@misc$SCdetMito_mito_ratio_info$input_scale, "percent")

  recorded_percent <- seu
  recorded_percent@misc$SCdetMito_mito_ratio_info <- list(
    mito_col = "mitoPercent",
    scale = "percent"
  )
  provenance_result <- ensure_mito_ratio(
    recorded_percent,
    mito_col = "mitoPercent",
    mito_scale = "auto",
    verbose = FALSE
  )
  expect_equal(provenance_result$mitoPercent, seu$mitoPercent / 100)

  expect_warning(
    ensure_mito_ratio(seu, mito_col = "mitoPercent", mito_scale = "auto", verbose = FALSE),
    "name suggests percentages"
  )

  seu$mitoRatio <- c(rep(0.10, 19), 1.10)
  expect_error(
    ensure_mito_ratio(seu, mito_col = "mitoRatio", mito_scale = "auto", verbose = FALSE),
    "ambiguous scale"
  )
})

test_that("numeric mitochondrial cutoffs have an explicit scale escape hatch", {
  expect_warning(
    expect_equal(normalize_mito_cutoff_value(1, scale = "auto"), 1),
    "100%"
  )
  expect_equal(normalize_mito_cutoff_value(1, scale = "percent"), 0.01)
  expect_error(normalize_mito_cutoff_value(2, scale = "fraction"), "between 0 and 1")
})

test_that("literature and quantile fallbacks are reported but never auto-applied", {
  seu <- make_threshold_safety_seurat()
  seu$mitoRatio <- seq(0.02, 0.08, length.out = ncol(seu))

  with_prior <- SCdetMito(
    seu,
    sample_col = "sample",
    species = "human",
    min_drop_cells = 1000,
    max_cut = 0.30,
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  )
  without_prior <- SCdetMito(
    seu,
    sample_col = "sample",
    min_drop_cells = 1000,
    max_cut = 0.30,
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  )

  expect_equal(with_prior$sample_cutoff_summary$recommended_cutoff, 0.10)
  expect_equal(with_prior$sample_cutoff_summary$recommendation_source, "literature_prior_fallback")
  expect_true(with_prior$sample_cutoff_summary$fallback_used)
  expect_false(with_prior$sample_cutoff_summary$recommended_auto_apply_eligible)
  expect_false(with_prior$recommended_auto_apply_eligible)
  expect_true(is.na(with_prior$data_driven_cutoff))

  expect_match(without_prior$sample_cutoff_summary$recommendation_source, "quantile_fallback")
  expect_true(without_prior$sample_cutoff_summary$fallback_used)
  expect_false(without_prior$sample_cutoff_summary$recommended_auto_apply_eligible)
})

test_that("upper search-boundary cutoffs require review", {
  seu <- make_threshold_safety_seurat()
  seu$mitoRatio <- c(rep(0.10, 20), rep(0.995, 10))

  result <- SCdetMito(
    seu,
    sample_col = "sample",
    min_cut = 0.01,
    max_cut = 1,
    bin_width = 0.01,
    min_drop_cells = 1,
    min_drop_fraction = 0,
    min_cells_after = 1,
    min_retention_after = 0,
    loss_test = "threshold_only",
    sample_cutoff_method = "first_significant_high",
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  )

  expect_equal(result$sample_cutoff_summary$selected_cutoff, 0.99)
  expect_true(result$sample_cutoff_summary$selected_upper_boundary_hit)
  expect_false(result$sample_cutoff_summary$selected_auto_apply_eligible)
  expect_equal(result$sample_cutoff_summary$recommendation_level, "review_required")
})

test_that("auto-apply guardrails flag evidence risk rather than imposing a universal low cap", {
  cautious <- build_recommendation_fields(
    first_significant_cutoff_high = 0.25,
    largest_drop_cutoff = 0.25,
    fallback_cutoff = NA_real_,
    fallback_used = FALSE,
    fallback_method = "none",
    fallback_quantile = 0.9,
    reference_cutoff = 0.10,
    retention_fraction_at_recommended = 0.80
  )
  extreme_reference <- build_recommendation_fields(
    first_significant_cutoff_high = 0.31,
    largest_drop_cutoff = 0.31,
    fallback_cutoff = NA_real_,
    fallback_used = FALSE,
    fallback_method = "none",
    fallback_quantile = 0.9,
    reference_cutoff = 0.10,
    retention_fraction_at_recommended = 0.80
  )
  exact_boundary <- build_recommendation_fields(
    first_significant_cutoff_high = 0.30,
    largest_drop_cutoff = 0.30,
    fallback_cutoff = NA_real_,
    fallback_used = FALSE,
    fallback_method = "none",
    fallback_quantile = 0.9,
    reference_cutoff = 0.10,
    retention_fraction_at_recommended = 0.80
  )
  low_retention <- build_recommendation_fields(
    first_significant_cutoff_high = 0.10,
    largest_drop_cutoff = 0.10,
    fallback_cutoff = NA_real_,
    fallback_used = FALSE,
    fallback_method = "none",
    fallback_quantile = 0.9,
    reference_cutoff = 0.10,
    retention_fraction_at_recommended = 0.29
  )
  below_auto_floor <- build_recommendation_fields(
    first_significant_cutoff_high = 0.10,
    largest_drop_cutoff = 0.10,
    fallback_cutoff = NA_real_,
    fallback_used = FALSE,
    fallback_method = "none",
    fallback_quantile = 0.9,
    reference_cutoff = 0.10,
    retention_fraction_at_recommended = 0.79
  )
  high_without_prior <- build_recommendation_fields(
    first_significant_cutoff_high = 0.50,
    largest_drop_cutoff = 0.50,
    fallback_cutoff = NA_real_,
    fallback_used = FALSE,
    fallback_method = "none",
    fallback_quantile = 0.9,
    reference_cutoff = NA_real_,
    retention_fraction_at_recommended = 0.80
  )

  expect_equal(cautious$recommendation_level, "cautious")
  expect_true(cautious$auto_apply_eligible)
  expect_equal(extreme_reference$recommendation_level, "review_required")
  expect_false(extreme_reference$auto_apply_eligible)
  expect_equal(exact_boundary$recommendation_status, "review_required")
  expect_false(exact_boundary$auto_apply_eligible)
  expect_false(low_retention$auto_apply_eligible)
  expect_equal(below_auto_floor$recommendation_level, "review_required")
  expect_false(below_auto_floor$auto_apply_eligible)
  expect_true(high_without_prior$auto_apply_eligible)
})

test_that("low-retention reference guidance yields an operating-domain alternative", {
  guarded <- build_recommendation_fields(
    first_significant_cutoff_high = 0.22,
    largest_drop_cutoff = 0.05,
    fallback_cutoff = NA_real_,
    fallback_used = FALSE,
    fallback_method = "none",
    fallback_quantile = 0.9,
    reference_cutoff = 0.05,
    retention_fraction_at_recommended = 0.82,
    reference_guided_cutoff = 0.05,
    reference_guided_retention = 0.10,
    first_significant_high_retention = 0.82,
    reference_retention = 0.10
  )

  expect_equal(guarded$recommended_cutoff, 0.22)
  expect_equal(guarded$recommended_method, "first_significant_high_retention_guard")
  expect_equal(
    guarded$recommendation_source,
    "reference_guided_overfilter_guard_first_significant_high"
  )
  expect_equal(guarded$recommendation_status, "review_required")
  expect_equal(guarded$recommendation_level, "review_required")
  expect_false(guarded$auto_apply_eligible)
  expect_match(guarded$recommended_reason, "at least 80%")
})

test_that("candidate and reference outside the operating domain yield no call", {
  no_call <- build_recommendation_fields(
    first_significant_cutoff_high = 0.22,
    largest_drop_cutoff = 0.05,
    fallback_cutoff = NA_real_,
    fallback_used = FALSE,
    fallback_method = "none",
    fallback_quantile = 0.9,
    reference_cutoff = 0.05,
    retention_fraction_at_recommended = 0.38,
    reference_guided_cutoff = 0.05,
    reference_guided_retention = 0.10,
    first_significant_high_retention = 0.38,
    reference_retention = 0.10
  )

  expect_true(is.na(no_call$recommended_cutoff))
  expect_equal(no_call$review_cutoff, 0.05)
  expect_equal(no_call$recommendation_status, "no_call")
  expect_equal(no_call$operating_domain_status, "outside_routine_domain")
  expect_false(no_call$auto_apply_eligible)
})

test_that("SCdetMito no-call is explicit and wrappers do not substitute it", {
  seu <- make_threshold_safety_seurat(n_cells = 30L)
  seu$mitoRatio <- c(rep(0.02, 6), rep(0.12, 24))

  detection <- SCdetMito(
    seu,
    sample_col = "sample",
    species = "human",
    min_cut = 0.01,
    max_cut = 0.20,
    bin_width = 0.01,
    min_drop_cells = 1,
    min_drop_fraction = 0,
    min_cells_after = 1,
    min_retention_after = 0,
    loss_test = "threshold_only",
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  )

  expect_equal(detection$sample_cutoff_summary$recommendation_status, "no_call")
  expect_true(is.na(detection$sample_cutoff_summary$recommended_cutoff))
  expect_true(is.finite(detection$sample_cutoff_summary$review_cutoff))
  expect_equal(detection$recommendation_status, "no_call")

  expect_error(
    SCQCone(
      seu,
      min_genes = 0,
      min_counts = 0,
      max_mito = "SCdetMito",
      species = "human",
      removeDouble = FALSE,
      plot = FALSE,
      table_out = FALSE
    ),
    "no actionable recommended cutoff"
  )
})

test_that("QC wrappers stop on review-only adaptive cutoffs unless explicitly overridden", {
  seu <- make_threshold_safety_seurat(n_cells = 10L)
  seu$mitoRatio <- rep(0.05, ncol(seu))

  expect_error(
    SCQCone(
      seu,
      min_genes = 0,
      min_counts = 0,
      max_mito = "SCdetMito",
      species = "human",
      removeDouble = FALSE,
      plot = FALSE,
      table_out = FALSE
    ),
    "marked review-only"
  )
  expect_error(
    SCQCmulti(
      seu,
      sample_col = "sample",
      mode = "all",
      min_genes = 0,
      min_counts = 0,
      max_mito = "SCdetMito",
      species = "human",
      removeDouble = FALSE,
      plot = FALSE,
      table_out = FALSE
    ),
    "marked review-only"
  )
  expect_warning(
    filtered <- SCQCone(
      seu,
      min_genes = 0,
      min_counts = 0,
      max_mito = "SCdetMito",
      species = "human",
      review_action = "warn_apply",
      removeDouble = FALSE,
      plot = FALSE,
      table_out = FALSE
    ),
    "marked review-only"
  )
  expect_s4_class(filtered, "Seurat")
  expect_false(filtered@misc$SCdetMito_QC$auto_apply_eligible)
})

test_that("benchmark method labels apply their selected cutoffs and exclude unsafe winners", {
  strategies <- normalize_benchmark_strategies(NULL, group_by = NULL)
  method_names <- c(
    "SCdetMito_largest_drop",
    "SCdetMito_first_significant_high",
    "SCdetMito_reference_guided"
  )
  method_strategies <- strategies[vapply(strategies, function(x) x$name %in% method_names, logical(1))]
  expect_true(all(!vapply(method_strategies, `[[`, logical(1), "use_recommended_cutoff")))

  score_df <- data.frame(
    strategy = c("SCdetMito_review_only", "fixed_10"),
    overall_score = c(1, 0.5),
    retention_rate = c(1, 0.9),
    median_mito_after = c(0.9, 0.05),
    sample_retention_cv = c(0, 0.1),
    group_retention_cv = c(NA_real_, NA_real_),
    auto_apply_eligible = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  recommendation <- recommend_benchmark_strategy(score_df)
  expect_true(is.na(recommendation$recommended_adaptive_strategy))
  expect_equal(recommendation$recommended_strategy, "fixed_10")
})
