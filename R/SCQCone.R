# SCdetMito
# Author: Silu Hu
# Contact: husilu0902@gmail.com
# Version: 1.4.1
# Last updated: 2026-05-21

#' SCQCone: perform QC for a single single-cell RNA-seq sample or group
#'
#' @description
#' `SCQCone()` performs quality control for a single Seurat object using
#' user-supplied thresholds or a mitochondrial cutoff inferred by `SCdetMito()`.
#'
#' @details
#' The function summarizes cell counts before and after QC, optionally exports
#' QC plots and tables, and can optionally run DoubletFinder-based doublet
#' filtering after the primary mitochondrial, gene-count, and UMI-count filters
#' are applied. Doublet removal is disabled by default because it is not part of
#' the core mitochondrial cutoff detector and requires the optional
#' `DoubletFinder` package.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data for a
#'   single sample.
#' @param feature_col,nFeature_RNA Metadata column storing detected feature
#'   counts. `feature_col` is preferred; `nFeature_RNA` is retained for
#'   backward compatibility.
#' @param count_col,nCount_RNA Metadata column storing UMI/count totals.
#'   `count_col` is preferred.
#' @param mito_col,mitoRatio Metadata column storing mitochondrial ratios.
#'   `mito_col` is preferred.
#' @param min_genes Minimum number of detected genes. Defaults to `200`.
#' @param max_genes Maximum number of detected genes. Defaults to `Inf`.
#' @param min_counts Minimum number of counts. Defaults to `500`.
#' @param max_counts Maximum number of counts. Defaults to `Inf`.
#' @param min_mito Minimum mitochondrial ratio. Defaults to `0`. Values between
#'   `1` and `100` are interpreted as percentages and converted to fractions.
#' @param max_mito Maximum mitochondrial ratio. Use `"SCdetMito"` to infer a
#'   cutoff automatically. Numeric values between `1` and `100` are interpreted
#'   as percentages and converted to fractions.
#' @param remove_doublets,removeDouble Whether to run optional
#'   DoubletFinder-based doublet removal after QC filtering. `remove_doublets`
#'   is preferred. Defaults to `FALSE`.
#' @param write_plots,plot Whether to export QC plots. `write_plots` is
#'   preferred. Defaults to `TRUE`.
#' @param write_tables,table_out Whether to export QC summary tables.
#'   `write_tables` is preferred. Defaults to `TRUE`.
#' @param output_dir Directory for exported files. Defaults to `"."`.
#' @param fail_action What to do if optional DoubletFinder removal fails.
#'   `"warn"` returns the mito-QC-filtered object with a warning; `"stop"`
#'   raises an error.
#' @param species,tissue Optional species and tissue names forwarded to
#'   [SCdetMito()] when `max_mito = "SCdetMito"` for reference-aware cutoff
#'   interpretation.
#' @param reference_cutoff,reference_table,reference_warning Optional reference
#'   cutoff controls forwarded to [SCdetMito()].
#' @param auto_add_mito Whether to calculate `mito_col` with
#'   [ensure_mito_ratio()] when the column is absent. Defaults to `TRUE`.
#' @param mito_features Optional complete mitochondrial feature vector forwarded
#'   to [ensure_mito_ratio()] and [SCdetMito()].
#' @param mito_pattern Optional mitochondrial feature regex forwarded to
#'   [ensure_mito_ratio()] and [SCdetMito()].
#' @param mito_assay Optional assay used for mitochondrial ratio calculation.
#' @param recompute_mito Whether to recompute `mito_col` even when it already
#'   exists. Defaults to `FALSE`.
#' @param ... Additional parameters, currently unused.
#'
#' @return A QC-filtered Seurat object.
#' @export
#'
#' @examples
#' # DO NOT RUN
#' # library(Seurat)
#' # counts <- matrix(rpois(2000, 5), nrow = 100)
#' # mt_genes <- c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2")
#' # # Representative toy mitochondrial genes only; not a complete gene set.
#' # rownames(counts) <- c(mt_genes, paste0("Gene", seq_len(96)))
#' # seu <- CreateSeuratObject(counts)
#' # seu <- ensure_mito_ratio(seu, species = "human")
#' # qc_seu <- SCQCone(seu, max_mito = 0.2, removeDouble = FALSE, plot = FALSE)
#' # qc_seu
#' # DO NOT RUN
SCQCone <- function(seurat_obj,
                    feature_col = NULL,
                    count_col = NULL,
                    mito_col = NULL,
                    nFeature_RNA = "nFeature_RNA",
                    nCount_RNA = "nCount_RNA",
                    mitoRatio = "mitoRatio",
                    min_genes = 200,
                    max_genes = Inf,
                    min_counts = 500,
                    max_counts = Inf,
                    min_mito = 0,
                    max_mito = "SCdetMito",
                    remove_doublets = NULL,
                    removeDouble = FALSE,
                    write_plots = NULL,
                    plot = TRUE,
                    write_tables = NULL,
                    table_out = TRUE,
                    output_dir = ".",
                    fail_action = c("warn", "stop"),
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
                    ...) {
  message("Performing quality control for single-cell RNA-seq data...")

  nFeature_RNA <- feature_col %||% nFeature_RNA
  nCount_RNA <- count_col %||% nCount_RNA
  mitoRatio <- mito_col %||% mitoRatio
  removeDouble <- remove_doublets %||% removeDouble
  plot <- write_plots %||% plot
  table_out <- write_tables %||% table_out
  fail_action <- match.arg(fail_action)

  output_dir <- ensure_output_dir(output_dir)
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
  qc_bounds <- validate_qc_bounds(
    min_genes = min_genes,
    max_genes = max_genes,
    min_counts = min_counts,
    max_counts = max_counts,
    min_mito = min_mito,
    max_mito = if (identical(max_mito, "SCdetMito")) NULL else max_mito
  )
  min_genes <- qc_bounds$min_genes
  max_genes <- qc_bounds$max_genes
  min_counts <- qc_bounds$min_counts
  max_counts <- qc_bounds$max_counts
  min_mito <- qc_bounds$min_mito

  inferred_mito_cutoff <- max_mito
  detection_details <- NULL
  if (identical(max_mito, "SCdetMito")) {
    temp_group_column <- ".scdetmito_temp_group"
    seurat_obj@meta.data[[temp_group_column]] <- rep("A", ncol(seurat_obj))
    detection_details <- SCdetMito(
      seurat_obj,
      mito_col = mitoRatio,
      sample_col = temp_group_column,
      write_tables = table_out,
      write_plots = plot,
      output_dir = output_dir,
      species = species,
      tissue = tissue,
      reference_cutoff = reference_cutoff,
      reference_table = reference_table,
      reference_warning = reference_warning,
      auto_add_mito = auto_add_mito,
      mito_features = mito_features,
      mito_pattern = mito_pattern,
      mito_assay = mito_assay,
      recompute_mito = recompute_mito,
      return_details = TRUE
    )
    inferred_mito_cutoff <- detection_details$cutoff
    seurat_obj@meta.data[[temp_group_column]] <- NULL
  } else {
    inferred_mito_cutoff <- qc_bounds$max_mito
  }

  reference_cell_count <- ncol(seurat_obj)
  before_summary <- build_qc_summary(
    seurat_obj,
    nFeature_RNA = nFeature_RNA,
    nCount_RNA = nCount_RNA,
    mitoRatio = mitoRatio
  )
  before_summary$Stage <- "before"
  before_summary$MitoCutoff <- inferred_mito_cutoff
  before_summary <- add_qc_outcome_metrics(before_summary, reference_cell_count)

  if (plot) {
    plot_obj <- seurat_obj
    plot_obj@meta.data$.scdetmito_temp_group <- rep("seurat_obj", ncol(plot_obj))
    SCQC_processedPlots(
      plot_obj,
      by = ".scdetmito_temp_group",
      flag = "Checked",
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      mitoRatio = mitoRatio,
      output_dir = output_dir
    )
  }

  message("QC Summary [before]:")
  print(before_summary)
  message("mitoRatio cutoff was set to: ", inferred_mito_cutoff)
  message("Other QC indicators:")
  message("min_genes: ", min_genes)
  message("max_genes: ", max_genes)
  message("min_counts: ", min_counts)
  message("max_counts: ", max_counts)

  message("Perform cell filtering ...")
  filtered_obj <- filter_cells_by_metrics(
    seurat_obj,
    nFeature_RNA = nFeature_RNA,
    nCount_RNA = nCount_RNA,
    mitoRatio = mitoRatio,
    min_genes = min_genes,
    max_genes = max_genes,
    min_counts = min_counts,
    max_counts = max_counts,
    min_mito = min_mito,
    max_mito = inferred_mito_cutoff
  )

  doublet_status <- list(
    requested = isTRUE(removeDouble),
    status = if (isTRUE(removeDouble)) "skipped" else "not_requested",
    error_message = NA_character_
  )
  if (isTRUE(removeDouble)) {
    message("Perform double-cell filtering ...")
    doublet_result <- run_doublet_filter_safely(
      filtered_obj,
      requested = TRUE,
      fail_action = fail_action
    )
    filtered_obj <- doublet_result$object
    doublet_status <- doublet_result$status
  }

  after_summary <- build_qc_summary(
    filtered_obj,
    nFeature_RNA = nFeature_RNA,
    nCount_RNA = nCount_RNA,
    mitoRatio = mitoRatio
  )
  after_summary$Stage <- "after"
  after_summary$MitoCutoff <- inferred_mito_cutoff
  after_summary <- add_qc_outcome_metrics(after_summary, reference_cell_count)

  if (plot) {
    plot_obj <- filtered_obj
    plot_obj@meta.data$.scdetmito_temp_group <- rep("seurat_obj", ncol(plot_obj))
    SCQC_processedPlots(
      plot_obj,
      by = ".scdetmito_temp_group",
      flag = "Filtered",
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      mitoRatio = mitoRatio,
      output_dir = output_dir
    )
  }

  message("QC Summary [after]:")
  print(after_summary)

  if (table_out) {
    qc_summary <- rbind(before_summary, after_summary)
    utils::write.csv(
      qc_summary,
      file = file.path(output_dir, "SCQCone_summary.csv"),
      row.names = FALSE
    )
  }

  fallback_events <- if (!is.null(detection_details$sample_cutoff_summary)) {
    detection_details$sample_cutoff_summary[
      detection_details$sample_cutoff_summary$fallback_used %in% TRUE,
      ,
      drop = FALSE
    ]
  } else {
    NULL
  }
  filtered_obj@misc$SCdetMito_QC <- build_qc_provenance(
    function_name = "SCQCone",
    input_cell_count = reference_cell_count,
    output_cell_count = ncol(filtered_obj),
    parameters = list(
      feature_col = nFeature_RNA,
      count_col = nCount_RNA,
      mito_col = mitoRatio,
      min_genes = min_genes,
      max_genes = max_genes,
      min_counts = min_counts,
      max_counts = max_counts,
      min_mito = min_mito,
      max_mito = max_mito,
      inferred_mito_cutoff = inferred_mito_cutoff,
      species = species,
      tissue = tissue,
      reference_cutoff = reference_cutoff,
      reference_warning = reference_warning,
      remove_doublets = removeDouble
    ),
    cutoff_plan = detection_details,
    cutoff_source = if (!is.null(detection_details$sample_cutoff_summary)) {
      paste(unique(detection_details$sample_cutoff_summary$cutoff_source), collapse = ";")
    } else {
      "user_defined"
    },
    cutoff_confidence = if (!is.null(detection_details$sample_cutoff_summary)) {
      paste(unique(detection_details$sample_cutoff_summary$cutoff_confidence), collapse = ";")
    } else {
      "user_defined"
    },
    loss_test = if (!is.null(detection_details)) detection_details$settings$loss_test %||% NA_character_ else NA_character_,
    sample_cutoff_method = if (!is.null(detection_details)) detection_details$settings$sample_cutoff_method %||% NA_character_ else NA_character_,
    fallback_events = fallback_events,
    doublet_filtering = doublet_status
  )

  filtered_obj
}
