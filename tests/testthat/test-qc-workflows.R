make_test_seurat <- function() {
  set.seed(42)
  counts <- matrix(
    rpois(200 * 80, lambda = 5),
    nrow = 200,
    ncol = 80
  )
  rownames(counts) <- c(paste0("MT-", 1:10), paste0("Gene", 11:200))
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))

  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  seu$sample <- rep(c("A", "B"), each = ncol(seu) / 2)
  seu$mitoRatio <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-") / 100
  seu
}

test_that("check_seu rescales percentage mitochondrial ratios", {
  seu <- make_test_seurat()
  seu$mitoPercent <- seu$mitoRatio * 100

  normalized <- check_seu(
    seu,
    check = "mitoPercent",
    normalize_fraction = TRUE,
    must_be_numeric = TRUE
  )

  expect_lte(max(normalized$mitoPercent, na.rm = TRUE), 1)
})

test_that("add_mitoRatio supports explicit livestock and bird species presets", {
  seu <- make_test_seurat()
  seu <- add_mitoRatio(seu, species = "pig", column = "mito_pig")
  seu <- add_mitoRatio(seu, species = "chicken", column = "mito_chicken")

  expect_true(all(c("mito_pig", "mito_chicken") %in% colnames(seu@meta.data)))
  expect_equal(seu$mito_pig, seu$mito_chicken)
})

test_that("add_mitoRatio detects pig-style bare mitochondrial gene symbols", {
  counts <- matrix(
    c(10, 5, 20, 15, 80, 90),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("ATP6", "ND1", "ACTB")
  colnames(counts) <- c("Cell1", "Cell2")

  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0))
  updated <- add_mitoRatio(seu, species = "pig", column = "mitoRatio")

  expected <- c((10 + 20) / (10 + 20 + 80), (5 + 15) / (5 + 15 + 90))
  expect_equal(unname(updated$mitoRatio), expected)
})

test_that("add_mitoRatio errors clearly when no mitochondrial genes are detected", {
  counts <- matrix(rpois(20 * 8, lambda = 5), nrow = 20, ncol = 8)
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0))

  expect_error(
    add_mitoRatio(seu, pattern = "^MT-", column = "mitoRatio"),
    "No mitochondrial features were detected"
  )
})

test_that("SCdetMito exports interval-loss tables and legacy aliases", {
  seu <- make_test_seurat()
  output_dir <- tempfile("scdetmito-det-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  cutoff <- SCdetMito(
    seu,
    by = "sample",
    plot = FALSE,
    table_out = TRUE,
    output_dir = output_dir
  )

  expect_type(cutoff, "double")
  expect_true(file.exists(file.path(output_dir, "SCdetMito_interval_loss_table.csv")))
  expect_true(file.exists(file.path(output_dir, "processed_temp_cell_num.csv")))
  expect_true(file.exists(file.path(output_dir, "mito_change_point_results.csv")))
})

test_that("SCdetMito default plotting works on demo object", {
  seu <- load_demo_seurat()
  output_dir <- tempfile("scdetmito-plot-test-")
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  det <- SCdetMito(
    seu,
    sample_col = "sample",
    mito_col = "mitoRatio",
    return_details = TRUE,
    write_plots = TRUE,
    output_dir = output_dir
  )

  expect_true(file.exists(file.path(output_dir, "SCdetMito_retained_cell_profiles.pdf")))
  expect_true(file.exists(file.path(output_dir, "SCdetMito_interval_loss_profiles.pdf")))
  expect_true(file.exists(file.path(output_dir, "SCdetMito_sample_cutoff_summary.pdf")))
  expect_true(is.data.frame(det$sample_cutoff_summary))
})

test_that("SCdetMito reference-guided plotting works on demo object", {
  seu <- load_demo_seurat()
  output_dir <- tempfile("scdetmito-ref-plot-test-")
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  det <- SCdetMito(
    seu,
    sample_col = "sample",
    mito_col = "mitoRatio",
    species = "human",
    tissue = "global",
    sample_cutoff_method = "reference_guided",
    return_details = TRUE,
    write_plots = TRUE,
    output_dir = output_dir
  )

  expect_true(file.exists(file.path(output_dir, "SCdetMito_retained_cell_profiles.pdf")))
  expect_true(file.exists(file.path(output_dir, "SCdetMito_interval_loss_profiles.pdf")))
  expect_true(file.exists(file.path(output_dir, "SCdetMito_sample_cutoff_summary.pdf")))
  expect_true("reference_cutoff" %in% names(det$sample_cutoff_summary))
})

test_that("SCdetMito supports alternative tests and adjustment methods", {
  seu <- make_test_seurat()

  result <- SCdetMito(
    seu,
    by = "sample",
    plot = FALSE,
    table_out = FALSE,
    loss_test = "empirical_tail",
    p_adjust_method = "none",
    sample_cutoff_method = "largest_drop",
    return_details = TRUE
  )

  expect_type(result$cutoff, "double")
  expect_equal(result$settings$loss_test, "empirical_tail")
  expect_equal(result$settings$p_adjust_method, "none")
  expect_equal(result$settings$detector_mode, "unadjusted_sensitivity_mode")
  expect_equal(result$settings$sample_cutoff_method, "largest_drop")
  expect_true(all(c("p_adjusted", "p_adjust_method") %in% colnames(result$change_points)))
  expect_identical(result$change_points, result$interval_loss_table)
  expect_true(all(c(
    "retained_cell_profile",
    "interval_loss_table",
    "significant_intervals",
    "sample_cutoff_summary",
    "settings"
  ) %in% names(result)))
  expect_true(all(c(
    "retained_cells",
    "interval_loss",
    "interval_loss_fraction",
    "retention_fraction",
    "passes_candidate_filters"
  ) %in% colnames(result$interval_loss_table)))
  expect_true("cutoff_confidence" %in% colnames(result$sample_cutoff_summary))
})

test_that("empirical-tail loss enrichment uses the same-sample empirical background", {
  result <- SCdetMito:::run_loss_enrichment_test(
    loss_value = 5,
    background_losses = c(1, 5, 7),
    method = "empirical_tail"
  )

  expect_equal(result$p_value, 0.75)
  expect_equal(result$test_method, "empirical_tail")
})

test_that("SCdetMito exposes selectable detection methods", {
  methods <- SCdetMito_methods()

  expect_true("loss_tests" %in% names(methods))
  expect_true(any(methods$loss_tests$method == "mad_zscore" & methods$loss_tests$default))
  expect_false(any(methods$loss_tests$method == "empirical_tail" & methods$loss_tests$default))
  expect_true("poisson_tail" %in% methods$loss_tests$method)
  expect_true("mad_zscore" %in% methods$loss_tests$method)
  expect_true("BH" %in% methods$p_adjust_methods$method)
  expect_true("reference_guided" %in% methods$sample_cutoff_methods$method)
  expect_true("first_significant_high" %in% methods$sample_cutoff_methods$method)
  expect_true("median_significant" %in% methods$sample_cutoff_methods$method)
  expect_true("min_significant" %in% methods$sample_cutoff_methods$method)
})

test_that("reference cutoff table is available", {
  ref <- SCdetMito_reference_cutoffs()

  expect_s3_class(ref, "data.frame")
  expect_true(all(c(
    "species",
    "tissue",
    "reference_cutoff",
    "reference_scope",
    "is_tissue_specific",
    "doi_or_url"
  ) %in% colnames(ref)))
  expect_true(any(ref$species == "human"))
  expect_true(any(ref$species == "mouse"))
  expect_true(all(ref$reference_cutoff > 0 & ref$reference_cutoff <= 1))
  expect_true(all(ref$doi_or_url == "https://doi.org/10.1093/bioinformatics/btaa751"))
  expect_false(any(ref$is_tissue_specific))
})

test_that("sample_cutoff_summary has unique column names", {
  seu <- load_demo_seurat()
  det <- SCdetMito(
    seu,
    sample_col = "sample",
    mito_col = "mitoRatio",
    return_details = TRUE,
    write_plots = FALSE
  )

  expect_equal(anyDuplicated(names(det$sample_cutoff_summary)), 0L)
  expect_true("n_cells" %in% names(det$sample_cutoff_summary))
})

test_that("detected_cutoff is backward-compatible alias of selected_cutoff", {
  seu <- load_demo_seurat()
  det <- SCdetMito(
    seu,
    sample_col = "sample",
    mito_col = "mitoRatio",
    species = "human",
    tissue = "global",
    sample_cutoff_method = "reference_guided",
    return_details = TRUE,
    write_plots = FALSE
  )
  summary <- det$sample_cutoff_summary

  expect_true("detected_cutoff" %in% names(summary))
  expect_true("selected_cutoff" %in% names(summary))
  expect_equal(summary$detected_cutoff, summary$selected_cutoff)
  expect_true(all(c(
    "has_reference_deviation",
    "has_low_retention_warning",
    "has_fallback_warning"
  ) %in% names(summary)))
})

test_that("empirical_tail with adjusted p-values warns about discreteness", {
  seu <- make_test_seurat()

  expect_warning(
    SCdetMito(
      seu,
      by = "sample",
      mitoRatio = "mitoRatio",
      loss_test = "empirical_tail",
      p_adj_method = "BH",
      plot = FALSE,
      table_out = FALSE,
      return_details = TRUE
    ),
    "rank-based and highly discrete"
  )
})

test_that("SCdetMito supports robust MAD z-score loss testing", {
  seu <- make_test_seurat()

  result <- suppressWarnings(SCdetMito(
    seu,
    by = "sample",
    mitoRatio = "mitoRatio",
    min_cut = 0.01,
    max_cut = 0.25,
    diff_num = 1,
    min_cells_after = 5,
    loss_test = "mad_zscore",
    p_adjust_method = "BH",
    alpha = 0.5,
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  ))

  expect_true(is.numeric(result$cutoff))
  expect_equal(result$settings$loss_test, "mad_zscore")
  expect_true("test_method" %in% colnames(result$change_points))
})

test_that("SCdetMito fallback none records not_detected sample cutoffs", {
  seu <- make_test_seurat()

  result <- suppressWarnings(SCdetMito(
    seu,
    by = "sample",
    mitoRatio = "mitoRatio",
    min_drop_cells = 10000,
    fallback_method = "none",
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  ))

  expect_true(all(is.na(result$sample_cutoff_summary$detected_cutoff)))
  expect_true(all(result$sample_cutoff_summary$cutoff_source == "not_detected"))
  expect_true(all(!result$sample_cutoff_summary$fallback_used))
})

test_that("largest_drop is the default sample-level cutoff selector", {
  counts <- matrix(rpois(50 * 100, lambda = 5), nrow = 50, ncol = 100)
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0))
  seu$sample <- "A"
  seu$mitoRatio <- c(rep(0.24, 20), rep(0.19, 5), rep(0.14, 15), rep(0.02, 60))

  result <- SCdetMito(
    seu,
    by = "sample",
    mitoRatio = "mitoRatio",
    bin_width = 0.05,
    min_cut = 0.05,
    max_cut = 0.25,
    min_drop_cells = 1,
    min_drop_fraction = 0,
    min_cells_after = 1,
    min_retention_after = 0,
    loss_test = "threshold_only",
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  )

  expect_equal(result$settings$sample_cutoff_method, "largest_drop")
  expect_equal(result$sample_cutoff_summary$detected_cutoff, 0.20)
  expect_equal(result$sample_cutoff_summary$selected_cutoff, result$sample_cutoff_summary$detected_cutoff)
  expect_true("largest_drop_cutoff" %in% colnames(result$sample_cutoff_summary))
  expect_equal(result$sample_cutoff_summary$cutoff_source, "significant_largest_drop")
})

test_that("reference-aware cutoff modes report first-boundary and selected cutoffs", {
  counts <- matrix(rpois(50 * 100, lambda = 5), nrow = 50, ncol = 100)
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0))
  seu$sample <- "A"
  seu$mitoRatio <- c(rep(0.24, 10), rep(0.07, 30), rep(0.02, 60))

  first <- SCdetMito(
    seu,
    by = "sample",
    mitoRatio = "mitoRatio",
    bin_width = 0.05,
    min_cut = 0.05,
    max_cut = 0.30,
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
  guided <- SCdetMito(
    seu,
    by = "sample",
    mitoRatio = "mitoRatio",
    bin_width = 0.05,
    min_cut = 0.05,
    max_cut = 0.30,
    min_drop_cells = 1,
    min_drop_fraction = 0,
    min_cells_after = 1,
    min_retention_after = 0,
    loss_test = "threshold_only",
    sample_cutoff_method = "reference_guided",
    species = "human",
    tissue = "global",
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  )

  expect_equal(first$sample_cutoff_summary$first_significant_cutoff_high, 0.20)
  expect_equal(first$sample_cutoff_summary$largest_drop_cutoff, 0.05)
  expect_equal(first$sample_cutoff_summary$detected_cutoff, first$sample_cutoff_summary$selected_cutoff)
  expect_equal(first$sample_cutoff_summary$cutoff_source, "first_significant_high")
  expect_equal(guided$sample_cutoff_summary$reference_cutoff, 0.10)
  expect_equal(guided$sample_cutoff_summary$selected_cutoff, 0.20)
  expect_true("warning_level" %in% colnames(guided$sample_cutoff_summary))
})

test_that("SCdetMito adjusts p-values within sample-level candidate intervals", {
  seu <- make_test_seurat()
  result <- SCdetMito(
    seu,
    by = "sample",
    mitoRatio = "mitoRatio",
    min_drop_cells = 1,
    min_drop_fraction = 0,
    min_cells_after = 1,
    min_retention_after = 0,
    loss_test = "mad_zscore",
    p_adj_method = "BH",
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  )
  candidates <- result$interval_loss_table[result$interval_loss_table$passes_candidate_filters, , drop = FALSE]

  expect_true(nrow(candidates) > 0)
  expect_true(all(is.finite(candidates$p_adj)))
  expect_true(all(candidates$p_adj_method == "BH"))
})

test_that("SCdetMito adjusts p-values across all interval-specific losses per sample", {
  seu <- make_test_seurat()
  result <- SCdetMito(
    seu,
    by = "sample",
    mitoRatio = "mitoRatio",
    min_drop_cells = 10000,
    loss_test = "mad_zscore",
    p_adj_method = "BH",
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  )
  interval_table <- result$interval_loss_table

  expect_true(all(is.finite(interval_table$p_adj)))
  expect_false(any(interval_table$passes_candidate_filters))
  expect_false(any(interval_table$is_significant))
})

test_that("SCdetMito supports p_adj_method none for statistical modes", {
  seu <- make_test_seurat()
  result <- SCdetMito(
    seu,
    by = "sample",
    mitoRatio = "mitoRatio",
    loss_test = "mad_zscore",
    p_adj_method = "none",
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  )
  interval_table <- result$interval_loss_table

  expect_equal(result$settings$p_adj_method, "none")
  expect_equal(result$settings$detector_mode, "unadjusted_sensitivity_mode")
  expect_equal(interval_table$p_adj, interval_table$p_value)
})

test_that("supported cutoff rule selects the largest majority-supported threshold", {
  support <- resolve_supported_cutoff(c(0.05, 0.08, 0.08, 0.20), min_support_fraction = 0.5)

  expect_equal(support$cutoff, 0.08)
  expect_equal(support$support_count, 3)
  expect_equal(support$required_support, 2)
})

test_that("SCdetMito final cutoff follows selected sample cutoffs in single-sample mode", {
  counts <- matrix(
    rpois(50 * 200, lambda = 5),
    nrow = 50,
    ncol = 200
  )
  rownames(counts) <- c(paste0("MT-", 1:5), paste0("Gene", 6:50))
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))

  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0))
  seu$sample <- "A"
  seu$mitoRatio <- rep(seq(0.05, 1.00, by = 0.05), each = 10)

  result <- SCdetMito(
    seu,
    by = "sample",
    mitoRatio = "mitoRatio",
    bin_width = 0.05,
    min_cut = 0.05,
    max_cut = 1.0,
    diff_num = 5,
    min_drop_fraction = 0,
    min_cells_after = 20,
    min_retention_after = 0,
    loss_test = "threshold_only",
    sample_cutoff_method = "median_significant",
    plot = FALSE,
    table_out = FALSE,
    return_details = TRUE
  )

  expect_equal(result$settings$sample_cutoff_method, "median_significant")
  expect_equal(result$cutoff, result$sample_cutoff_summary$detected_cutoff[[1]])
  expect_true("supports_final_cutoff" %in% colnames(result$sample_cutoff_summary))
})

test_that("SCQCone and SCQCmulti support custom metric columns", {
  seu <- make_test_seurat()
  seu$genes_detected <- seu$nFeature_RNA
  seu$umi_detected <- seu$nCount_RNA
  seu$mito_fraction <- seu$mitoRatio

  output_dir <- tempfile("scdetmito-qc-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  qc_one <- SCQCone(
    seu,
    nFeature_RNA = "genes_detected",
    nCount_RNA = "umi_detected",
    mitoRatio = "mito_fraction",
    min_genes = 0,
    min_counts = 0,
    max_mito = 0.2,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = TRUE,
    output_dir = output_dir
  )
  qc_all <- SCQCmulti(
    seu,
    by = "sample",
    mode = "all",
    nFeature_RNA = "genes_detected",
    nCount_RNA = "umi_detected",
    mitoRatio = "mito_fraction",
    min_genes = 0,
    min_counts = 0,
    max_mito = 0.2,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = TRUE,
    output_dir = output_dir
  )
  qc_split <- SCQCmulti(
    seu,
    by = "sample",
    mode = "split",
    nFeature_RNA = "genes_detected",
    nCount_RNA = "umi_detected",
    mitoRatio = "mito_fraction",
    min_genes = 0,
    min_counts = 0,
    max_mito = 0.2,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = TRUE,
    output_dir = output_dir
  )
  qc_adaptive <- SCQCmulti(
    seu,
    by = "sample",
    mode = "all",
    max_mito = "SCdetMito",
    scdet_options = list(
      loss_test = "threshold_only",
      p_adjust_method = "bonferroni",
      sample_cutoff_method = "median_significant"
    ),
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = TRUE,
    output_dir = output_dir
  )

  expect_s4_class(qc_one, "Seurat")
  expect_s4_class(qc_all, "Seurat")
  expect_s4_class(qc_split, "Seurat")
  expect_s4_class(qc_adaptive, "Seurat")
  expect_equal(ncol(qc_all), ncol(qc_split))
  expect_true(file.exists(file.path(output_dir, "SCQCone_summary.csv")))
  expect_true(file.exists(file.path(output_dir, "SCQCmulti_all_summary.csv")))
  expect_true(file.exists(file.path(output_dir, "SCQCmulti_split_summary.csv")))
  expect_equal(qc_adaptive@misc$SCdetMito_QC$detection$settings$p_adjust_method, "bonferroni")
  expect_equal(qc_one@misc$SCdetMito_QC$package, "SCdetMito")
  expect_true(is.list(qc_all@misc$SCdetMito_QC$provenance))
})

test_that("SCQCmulti adaptive strategies support consensus, strictest, and groupwise plans", {
  seu <- make_test_seurat()
  seu$group <- ifelse(seu$sample == "A", "control", "case")

  consensus <- SCQCmulti(
    seu,
    sample_col = "sample",
    group_col = "group",
    mode = "all",
    cutoff_strategy = "consensus",
    max_mito = "SCdetMito",
    scdet_options = list(loss_test = "threshold_only", min_drop_cells = 1),
    min_genes = 0,
    min_counts = 0,
    remove_doublets = FALSE,
    write_plots = FALSE,
    write_tables = FALSE
  )
  strictest <- SCQCmulti(
    seu,
    sample_col = "sample",
    group_col = "group",
    mode = "all",
    cutoff_strategy = "strictest",
    max_mito = "SCdetMito",
    scdet_options = list(loss_test = "threshold_only", min_drop_cells = 1),
    min_genes = 0,
    min_counts = 0,
    remove_doublets = FALSE,
    write_plots = FALSE,
    write_tables = FALSE
  )
  groupwise <- suppressWarnings(SCQCmulti(
    seu,
    sample_col = "sample",
    group_col = "group",
    mode = "all",
    cutoff_strategy = "groupwise",
    max_mito = "SCdetMito",
    scdet_options = list(loss_test = "threshold_only", min_drop_cells = 1),
    min_genes = 0,
    min_counts = 0,
    remove_doublets = FALSE,
    write_plots = FALSE,
    write_tables = FALSE
  ))

  expect_s4_class(consensus, "Seurat")
  expect_s4_class(strictest, "Seurat")
  expect_s4_class(groupwise, "Seurat")
  expect_equal(consensus@misc$SCdetMito_QC$strategy, "consensus")
  expect_equal(strictest@misc$SCdetMito_QC$strategy, "strictest")
  expect_equal(groupwise@misc$SCdetMito_QC$strategy, "groupwise")
})

test_that("SCQCmulti warns for groupwise cutoffs with single-sample groups", {
  seu <- make_test_seurat()
  seu$group <- ifelse(seu$sample == "A", "control", "case")

  expect_warning(
    SCQCmulti(
      seu,
      sample_col = "sample",
      group_col = "group",
      mode = "all",
      cutoff_strategy = "groupwise",
      max_mito = 0.2,
      min_genes = 0,
      min_counts = 0,
      remove_doublets = FALSE,
      write_plots = FALSE,
      write_tables = FALSE
    ),
    "fewer than two samples"
  )
})

test_that("SCdetMito_sensitivity returns sample-level cutoff stability rows", {
  seu <- make_test_seurat()

  sensitivity <- SCdetMito_sensitivity(
    seu,
    sample_col = "sample",
    mito_col = "mitoRatio",
    bin_width = c(0.02),
    alpha = c(0.05),
    min_drop_fraction = c(0),
    loss_test = c("threshold_only"),
    sample_cutoff_method = c("largest_drop", "max_significant"),
    min_drop_cells = 1,
    min_cells_after = 1,
    min_retention_after = 0,
    write_tables = FALSE,
    write_plots = FALSE
  )

  expect_s3_class(sensitivity, "data.frame")
  expect_true(all(c(
    "sample",
    "parameter_id",
    "detected_cutoff",
    "cutoff_source",
    "retained_cells_at_cutoff",
    "retention_fraction_at_cutoff"
  ) %in% colnames(sensitivity)))
  expect_equal(length(unique(sensitivity$sample_cutoff_method)), 2)

  sensitivity_details <- SCdetMito_sensitivity(
    seu,
    sample_col = "sample",
    mito_col = "mitoRatio",
    bin_width = c(0.02),
    alpha = c(0.05),
    min_drop_fraction = c(0),
    loss_test = c("threshold_only"),
    sample_cutoff_method = c("largest_drop"),
    min_drop_cells = 1,
    min_cells_after = 1,
    min_retention_after = 0,
    return_details = TRUE
  )
  expect_true(all(c(
    "cutoff_median",
    "cutoff_iqr",
    "cutoff_min",
    "cutoff_max",
    "fallback_fraction",
    "not_detected_fraction"
  ) %in% colnames(sensitivity_details$sensitivity_summary)))
})

test_that("QC functions normalize percentage-style mitochondrial cutoffs", {
  seu <- make_test_seurat()

  expect_warning(
    qc_percent <- SCQCone(
      seu,
      min_genes = 0,
      min_counts = 0,
      max_mito = 20,
      removeDouble = FALSE,
      plot = FALSE,
      table_out = FALSE
    ),
    "converted to a fraction"
  )
  qc_fraction <- SCQCone(
    seu,
    min_genes = 0,
    min_counts = 0,
    max_mito = 0.2,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = FALSE
  )

  expect_equal(ncol(qc_percent), ncol(qc_fraction))
  expect_error(
    SCQCone(
      seu,
      min_genes = 0,
      min_counts = 0,
      min_mito = 0.4,
      max_mito = 0.2,
      removeDouble = FALSE,
      plot = FALSE,
      table_out = FALSE
    ),
    "min_mito"
  )
})
