test_that("SCQCbenchmark produces comparative outputs on demo data", {
  demo_path <- system.file("extdata", "demo_multigroup_seurat.rds", package = "SCdetMito")
  seu <- readRDS(demo_path)
  output_dir <- tempfile("scdetmito-benchmark-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  strategies <- list(
    list(name = "fixed_0.10", mode = "all", cutoff_strategy = "consensus", max_mito = 0.10),
    list(
      name = "scdet_consensus",
      mode = "all",
      cutoff_strategy = "consensus",
      max_mito = "SCdetMito",
      scdet_options = list(loss_test = "empirical_tail")
    ),
    list(
      name = "scdet_groupwise",
      mode = "all",
      cutoff_strategy = "groupwise",
      max_mito = "SCdetMito",
      scdet_options = list(sample_cutoff_method = "largest_drop")
    )
  )

  benchmark <- suppressWarnings(SCQCbenchmark(
    seurat_obj = seu,
    sample_by = "sample",
    group_by = "group",
    strategies = strategies,
    scdet_options = list(p_adjust_method = "holm"),
    score_weights = c(
      retention = 0.30,
      mito = 0.40,
      sample_balance = 0.20,
      group_balance = 0.10
    ),
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    run_downstream = TRUE,
    plot = TRUE,
    save_objects = FALSE,
    output_dir = output_dir
  ))

  expect_true(file.exists(file.path(output_dir, "benchmark_summary.csv")))
  expect_true(file.exists(file.path(output_dir, "benchmark_sample_retention.csv")))
  expect_true(file.exists(file.path(output_dir, "benchmark_group_retention.csv")))
  expect_true(file.exists(file.path(output_dir, "benchmark_report.md")))
  expect_true(file.exists(file.path(output_dir, "benchmark_scores.csv")))
  expect_true(file.exists(file.path(output_dir, "benchmark_cells_retained.pdf")))
  expect_true(file.exists(file.path(output_dir, "benchmark_sample_retention.pdf")))
  expect_true(all(c("strategy", "retention_rate", "applied_cutoff_min") %in% colnames(benchmark$summary)))
  expect_true(all(c("strategy", "overall_score", "retention_score") %in% colnames(benchmark$scores)))
  expect_true("adaptive_strategy" %in% colnames(benchmark$recommendation))
  expect_true("recommended_adaptive_strategy" %in% colnames(benchmark$recommendation))
  expect_true("recommended_strategy" %in% colnames(benchmark$recommendation))
  expect_true("recommended_adaptive_strategy" %in% colnames(benchmark$recommended_strategies))
  expect_equal(
    benchmark$strategy_objects$scdet_consensus@misc$SCdetMito_QC$detection$settings$loss_test,
    "empirical_tail"
  )
  expect_equal(
    benchmark$strategy_objects$scdet_groupwise@misc$SCdetMito_QC$detection$settings$p_adjust_method,
    "holm"
  )

  summary_df <- benchmark$summary
  consensus_cells <- summary_df$total_cells_after[summary_df$strategy == "scdet_consensus"]
  groupwise_cells <- summary_df$total_cells_after[summary_df$strategy == "scdet_groupwise"]

  expect_true(length(consensus_cells) == 1)
  expect_true(length(groupwise_cells) == 1)
  expect_gte(groupwise_cells, consensus_cells)
  expect_true(benchmark$recommendation$adaptive_strategy %in% c("scdet_consensus", "scdet_groupwise"))
  expect_equal(
    benchmark$recommendation$adaptive_strategy,
    benchmark$recommendation$recommended_adaptive_strategy
  )
})

test_that("SCQCbenchmark handles single-sample adaptive benchmarks", {
  demo_path <- system.file("extdata", "demo_multigroup_seurat.rds", package = "SCdetMito")
  seu <- readRDS(demo_path)
  cells_to_keep <- rownames(subset(seu@meta.data, sample == "S1"))
  seu <- subset(seu, cells = cells_to_keep)
  output_dir <- tempfile("scdetmito-single-benchmark-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  benchmark <- suppressWarnings(SCQCbenchmark(
    seurat_obj = seu,
    sample_by = "sample",
    group_by = NULL,
    strategies = list(
      list(name = "fixed_0.10", mode = "all", cutoff_strategy = "consensus", max_mito = 0.10),
      list(name = "scdet_consensus", mode = "all", cutoff_strategy = "consensus", max_mito = "SCdetMito")
    ),
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    run_downstream = FALSE,
    plot = FALSE,
    output_dir = output_dir
  ))

  expect_true(all(is.finite(benchmark$summary$sample_retention_cv)))
  expect_true(all(is.finite(benchmark$scores$overall_score)))
  expect_equal(benchmark$recommendation$adaptive_strategy, "scdet_consensus")
})
