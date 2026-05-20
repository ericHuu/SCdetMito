test_that("demo Seurat RDS is available and loadable", {
  demo_path <- system.file("extdata", "demo_multigroup_seurat.rds", package = "SCdetMito")

  expect_true(nzchar(demo_path))

  seu <- readRDS(demo_path)
  expect_s4_class(seu, "Seurat")
  expect_equal(length(unique(seu$sample)), 6)
  expect_equal(length(unique(seu$group)), 3)
  expect_true("mitoRatio" %in% colnames(seu@meta.data))
})

test_that("load_demo_seurat returns the bundled Seurat object", {
  seu <- load_demo_seurat()

  expect_s4_class(seu, "Seurat")
  expect_equal(length(unique(seu$sample)), 6)
  expect_equal(length(unique(seu$group)), 3)
})

test_that("add_mitoRatio computes mitochondrial fractions", {
  demo_path <- system.file("extdata", "demo_multigroup_seurat.rds", package = "SCdetMito")
  seu <- readRDS(demo_path)
  seu$mitoRatio <- NULL

  updated <- add_mitoRatio(seu, column = "mitoRatio")

  expect_true("mitoRatio" %in% colnames(updated@meta.data))
  expect_lte(max(updated$mitoRatio, na.rm = TRUE), 1)
  expect_gte(min(updated$mitoRatio, na.rm = TRUE), 0)
})

test_that("add_mitoRatio auto-detects lowercase mouse-style mitochondrial genes", {
  counts <- matrix(
    c(10, 5, 20, 15, 80, 90),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("mt-Atp6", "mt-Nd1", "Actb")
  colnames(counts) <- c("Cell1", "Cell2")

  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0))
  updated <- add_mitoRatio(seu, column = "mitoRatio")

  expected <- c((10 + 20) / (10 + 20 + 80), (5 + 15) / (5 + 15 + 90))
  expect_equal(unname(updated$mitoRatio), expected)
})

test_that("add_mitoRatio fails clearly when no mitochondrial features are present", {
  counts <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    byrow = TRUE
  )
  rownames(counts) <- c("Actb", "Gapdh")
  colnames(counts) <- c("Cell1", "Cell2")

  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0))

  expect_error(
    add_mitoRatio(seu),
    "No mitochondrial features were detected"
  )
})

test_that("group-aware SCQCmulti strategies run on bundled demo data", {
  demo_path <- system.file("extdata", "demo_multigroup_seurat.rds", package = "SCdetMito")
  seu <- readRDS(demo_path)
  output_dir <- tempfile("scdetmito-demo-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  qc_consensus <- suppressWarnings(SCQCmulti(
    seu,
    by = "sample",
    group_by = "group",
    mode = "all",
    cutoff_strategy = "consensus",
    max_mito = "SCdetMito",
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = TRUE,
    output_dir = output_dir
  ))
  qc_groupwise <- suppressWarnings(SCQCmulti(
    seu,
    by = "sample",
    group_by = "group",
    mode = "all",
    cutoff_strategy = "groupwise",
    max_mito = "SCdetMito",
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = TRUE,
    output_dir = output_dir
  ))
  qc_split <- suppressWarnings(SCQCmulti(
    seu,
    by = "sample",
    group_by = "group",
    mode = "split",
    cutoff_strategy = "groupwise",
    max_mito = "SCdetMito",
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    plot = FALSE,
    table_out = TRUE,
    output_dir = output_dir
  ))

  plan_table <- utils::read.csv(file.path(output_dir, "SCQCmulti_cutoff_plan.csv"))
  group_table <- utils::read.csv(file.path(output_dir, "SCQCmulti_group_cutoffs.csv"))

  expect_s4_class(qc_consensus, "Seurat")
  expect_s4_class(qc_groupwise, "Seurat")
  expect_s4_class(qc_split, "Seurat")
  expect_equal(length(unique(qc_consensus@misc$SCdetMito_QC$sample_plan$applied_cutoff)), 1)
  expect_equal(qc_consensus@misc$SCdetMito_QC$group_plan$group_id, "all_samples")
  expect_equal(qc_consensus@misc$SCdetMito_QC$sample_plan$aggregation_level, rep("global", 6))
  expect_gt(length(unique(qc_groupwise@misc$SCdetMito_QC$sample_plan$applied_cutoff)), 1)
  expect_gte(ncol(qc_groupwise), ncol(qc_consensus))
  expect_equal(ncol(qc_groupwise), ncol(qc_split))
  expect_true(all(c("sample_id", "group_id", "applied_cutoff", "supports_final_cutoff") %in% colnames(plan_table)))
  expect_true(all(c("group_id", "group_cutoff", "support_count", "support_fraction") %in% colnames(group_table)))
  expect_true(file.exists(file.path(output_dir, "SCQCmulti_group_cutoffs.csv")))
})
