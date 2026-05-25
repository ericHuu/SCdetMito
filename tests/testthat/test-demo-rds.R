test_that("load_demo_pbmc returns the bundled multi-donor PBMC object", {
  demo_path <- system.file("extdata", "demo_pbmc_multidonor_seurat.rds", package = "SCdetMito")

  expect_true(nzchar(demo_path))

  seu <- load_demo_pbmc()
  expect_s4_class(seu, "Seurat")
  expect_equal(ncol(seu), 10000)
  expect_equal(length(unique(seu$sample)), 4)
  expect_equal(length(unique(seu$group)), 2)
  expect_true(all(c(
    "sample", "donor", "sex", "age_group", "group", "condition",
    "tissue", "species", "dataset", "data_source", "demo_note", "mitoRatio"
  ) %in% colnames(seu@meta.data)))
  expect_true(is.numeric(seu$mitoRatio))
  expect_lte(max(seu$mitoRatio, na.rm = TRUE), 1)
  expect_gte(min(seu$mitoRatio, na.rm = TRUE), 0)
})

test_that("load_demo_seurat is a backward-compatible demo alias", {
  expect_s4_class(load_demo_seurat(), "Seurat")
  expect_equal(ncol(load_demo_seurat()), ncol(load_demo_pbmc()))
})

test_that("bundled PBMC demo retains raw counts and mitochondrial genes", {
  seu <- load_demo_pbmc()
  assay <- Seurat::DefaultAssay(seu)
  mt_genes <- grep("^MT-", rownames(seu), value = TRUE)
  layers_fun <- get("Layers", envir = asNamespace("Seurat"))
  layer_data_fun <- get("LayerData", envir = asNamespace("Seurat"))
  count_layers <- grep(
    "^counts",
    layers_fun(seu[[assay]]),
    value = TRUE
  )

  expect_gt(length(mt_genes), 0)
  expect_gt(length(count_layers), 0)
  mt_count_sum <- sum(vapply(
    count_layers,
    function(layer) {
      layer_counts <- layer_data_fun(seu[[assay]], layer = layer)
      sum(layer_counts[intersect(mt_genes, rownames(layer_counts)), , drop = FALSE])
    },
    numeric(1)
  ))
  expect_gt(mt_count_sum, 0)
})

test_that("bundled PBMC demo mitoRatio recomputes consistently", {
  seu <- load_demo_pbmc()
  existing <- seu$mitoRatio

  seu2 <- load_demo_pbmc()
  seu2$mitoRatio <- NULL
  seu2 <- ensure_mito_ratio(seu2, mito_col = "mitoRatio", species = "human", verbose = FALSE)

  expect_equal(existing, seu2$mitoRatio, tolerance = 1e-8)
})

test_that("SCdetMito runs on bundled PBMC demo with reference guidance", {
  seu <- load_demo_pbmc()

  det <- SCdetMito(
    seu,
    sample_col = "sample",
    mito_col = "mitoRatio",
    species = "human",
    tissue = "PBMC",
    sample_cutoff_method = "reference_guided",
    return_details = TRUE,
    write_plots = FALSE,
    write_tables = FALSE
  )

  expect_true(is.data.frame(det$sample_cutoff_summary))
  expect_equal(anyDuplicated(names(det$sample_cutoff_summary)), 0L)
  expect_true(all(c("selected_cutoff", "detected_cutoff") %in% names(det$sample_cutoff_summary)))
  expect_equal(det$sample_cutoff_summary$selected_cutoff, det$sample_cutoff_summary$detected_cutoff)
})

test_that("SCdetMito cutoffs are consistent after demo mitoRatio recomputation", {
  existing <- SCdetMito(
    load_demo_pbmc(),
    sample_col = "sample",
    mito_col = "mitoRatio",
    species = "human",
    tissue = "PBMC",
    sample_cutoff_method = "reference_guided",
    return_details = TRUE,
    write_plots = FALSE,
    write_tables = FALSE
  )

  recomputed_obj <- load_demo_pbmc()
  recomputed_obj$mitoRatio <- NULL
  recomputed <- SCdetMito(
    recomputed_obj,
    sample_col = "sample",
    mito_col = "mitoRatio",
    species = "human",
    tissue = "PBMC",
    sample_cutoff_method = "reference_guided",
    auto_add_mito = TRUE,
    return_details = TRUE,
    write_plots = FALSE,
    write_tables = FALSE
  )

  expect_equal(
    existing$sample_cutoff_summary$selected_cutoff,
    recomputed$sample_cutoff_summary$selected_cutoff,
    tolerance = 1e-8
  )
})

test_that("SCQCmulti runs groupwise on bundled PBMC demo", {
  seu <- load_demo_pbmc()

  qc <- suppressWarnings(SCQCmulti(
    seu,
    sample_col = "sample",
    group_col = "group",
    mito_col = "mitoRatio",
    cutoff_strategy = "groupwise",
    scdet_options = list(sample_cutoff_method = "reference_guided"),
    species = "human",
    tissue = "PBMC",
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    write_plots = FALSE,
    write_tables = FALSE
  ))

  expect_s4_class(qc, "Seurat")
  expect_true(!is.null(qc@misc$SCdetMito_QC$sample_plan))
})

test_that("SCQCmulti groupwise provenance reports final strategy-level cutoffs", {
  seu <- load_demo_pbmc()

  qc_obj <- suppressWarnings(SCQCmulti(
    seu,
    sample_col = "sample",
    group_col = "group",
    mito_col = "mitoRatio",
    cutoff_strategy = "groupwise",
    scdet_options = list(sample_cutoff_method = "reference_guided"),
    species = "human",
    tissue = "PBMC",
    min_genes = 0,
    min_counts = 0,
    removeDouble = FALSE,
    write_plots = FALSE,
    write_tables = FALSE
  ))

  qc <- qc_obj@misc$SCdetMito_QC
  expect_true(is.data.frame(qc$sample_cutoff_summary))
  expect_true(is.data.frame(qc$group_cutoff_summary))
  expect_equal(qc$cutoff_applied_source, "recommended_cutoff")
  expect_equal(qc$cutoff_applied_strategy, "groupwise")
  expect_equal(qc$applied_cutoff_level, "group")
  expect_equal(qc$final_recommended_cutoff, qc$cutoff_applied)
  expect_equal(qc$recommended_cutoff, qc$cutoff_applied)
  expect_true("recommended_cutoff" %in% colnames(qc$sample_cutoff_summary))
  expect_true("group_cutoff" %in% colnames(qc$group_cutoff_summary))
  expect_equal(qc$cutoff_applied, max(qc$group_cutoff_summary$group_cutoff, na.rm = TRUE))
})
