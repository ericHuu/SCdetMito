# SCdetMito
# Author: Silu Hu
# Contact: husilu0902@gmail.com
# Version: 1.4.3
# Last updated: 2026-05-23

#' SCQCmulti: perform QC for multi-sample single-cell RNA-seq data
#'
#' @description
#' `SCQCmulti()` performs quality control for integrated multi-sample scRNA-seq
#' datasets. It supports both sample-wise QC and full-object QC with
#' sample-level or group-level mitochondrial cutoffs.
#'
#' @details
#' The function is designed for studies in which several samples may belong to
#' one biological group. When `max_mito = "SCdetMito"`, sample-level cutoffs
#' are first estimated and then aggregated according to the selected
#' `cutoff_strategy`. `SCQCmulti()` stores both sample-level cutoff evidence
#' and strategy-level applied cutoffs. The field `cutoff_applied` records the
#' actual filtering cutoff, while `sample_cutoff_summary` and
#' `group_cutoff_summary` retain the intermediate evidence used to derive it.
#' In groupwise mode, the top-level `recommended_cutoff` represents the final
#' strategy-level recommended cutoff, not an arbitrary sample-level value.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param sample_col,by Metadata column identifying samples. For multi-group
#'   studies this should be the sample column. `sample_col` is preferred; `by`
#'   is retained for backward compatibility.
#' @param group_col,group_by Optional metadata column defining higher-level
#'   biological groups that contain multiple samples. `group_col` is preferred.
#' @param mode QC strategy. `"split"` runs QC per sample and merges the
#'   results; `"all"` runs QC on the full object using cell-level cutoffs
#'   derived from the selected strategy.
#' @param cutoff_strategy Strategy used when `max_mito = "SCdetMito"`.
#'   `"consensus"` derives sample-level cutoffs, selects one all-sample
#'   sample-supported global cutoff for the integrated object, and applies it
#'   across cells even when `group_col`/`group_by` is provided.
#'   `"strictest"` applies the minimum sample-level cutoff across all samples.
#'   `"groupwise"` is the only strategy that applies group-specific cutoffs, or
#'   sample-specific cutoffs when `group_by = NULL`; the reported final cutoff
#'   is the largest group-level supported cutoff.
#' @param cutoff_quantile Quantile used to summarize sample-level cutoffs within
#'   a group for legacy workflows. Defaults to `0.25`.
#' @param cutoff_support_fraction Minimum fraction of samples required to support
#'   a global or group-level recommended cutoff. A sample supports a cutoff when
#'   its sample-level detected cutoff is at least that value. Defaults to `0.5`.
#' @param scdet_options Optional list of additional arguments forwarded to
#'   [SCdetMito()] when `max_mito = "SCdetMito"`. This can be used to select
#'   alternative loss tests, multiple-testing correction methods, and
#'   sample-level cutoff rules.
#' @param feature_col,nFeature_RNA Metadata column storing detected feature
#'   counts. `feature_col` is preferred.
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
#'   [SCdetMito()] when adaptive mitochondrial cutoffs are requested.
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
#' @param use_recommended_cutoff When `max_mito = "SCdetMito"`, use
#'   `recommended_cutoff` from [SCdetMito()] as the sample-level input for
#'   consensus, strictest, and groupwise filtering. Defaults to `TRUE`. Set to
#'   `FALSE` to apply the user-selected `selected_cutoff`.
#' @param ... Additional parameters, currently unused.
#'
#' @return A QC-filtered Seurat object.
#' @export
#'
#' @examples
#' # DO NOT RUN
#' # library(Seurat)
#' # counts <- matrix(rpois(4000, 5), nrow = 100)
#' # mt_genes <- c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2")
#' # # Representative toy mitochondrial genes only; not a complete gene set.
#' # rownames(counts) <- c(mt_genes, paste0("Gene", seq_len(96)))
#' # seu <- CreateSeuratObject(counts)
#' # seu$sample <- rep(c("A", "B"), each = ncol(seu) / 2)
#' # seu <- ensure_mito_ratio(seu, species = "human")
#' # qc_seu <- SCQCmulti(seu, by = "sample", mode = "split",
#' #   max_mito = 0.2, removeDouble = FALSE, plot = FALSE
#' # )
#' # qc_seu
#' # DO NOT RUN
SCQCmulti <- function(seurat_obj,
                      sample_col = NULL,
                      group_col = NULL,
                      by = NULL,
                      group_by = NULL,
                      mode = c("split", "all"),
                      cutoff_strategy = c("consensus", "strictest", "groupwise"),
                      cutoff_quantile = 0.25,
                      cutoff_support_fraction = 0.5,
                      scdet_options = list(),
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
                      use_recommended_cutoff = TRUE,
                      ...) {
  message("Performing quality control for single-cell RNA-seq data...")

  by <- sample_col %||% by
  group_by <- group_col %||% group_by
  nFeature_RNA <- feature_col %||% nFeature_RNA
  nCount_RNA <- count_col %||% nCount_RNA
  mitoRatio <- mito_col %||% mitoRatio
  removeDouble <- remove_doublets %||% removeDouble
  plot <- write_plots %||% plot
  table_out <- write_tables %||% table_out
  fail_action <- match.arg(fail_action)
  if (is.null(by) || identical(by, "")) {
    stop("A sample column must be supplied through 'sample_col' or legacy 'by'.", call. = FALSE)
  }

  mode <- match.arg(mode)
  cutoff_strategy <- match.arg(cutoff_strategy)
  output_dir <- ensure_output_dir(output_dir)
  group_by <- normalize_group_by(group_by)

  seurat_obj <- check_seu(seurat_obj, by)
  if (!is.null(group_by)) {
    seurat_obj <- check_seu(seurat_obj, group_by)
    if (identical(cutoff_strategy, "groupwise")) {
      sample_group_pairs <- unique(seurat_obj@meta.data[, c(by, group_by), drop = FALSE])
      group_sample_counts <- table(as.character(sample_group_pairs[[group_by]]))
      small_groups <- names(group_sample_counts)[group_sample_counts < 2L]
      if (length(small_groups)) {
        warning(
          "Some groups contain fewer than two samples. Groupwise sample-supported cutoffs are less stable for single-sample groups.",
          call. = FALSE
        )
      }
    }
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
  if (anyNA(seurat_obj@meta.data[[by]])) {
    stop("The grouping column specified by 'by' cannot contain missing values.", call. = FALSE)
  }
  if (!is.null(group_by) && anyNA(seurat_obj@meta.data[[group_by]])) {
    stop("The grouping column specified by 'group_by' cannot contain missing values.", call. = FALSE)
  }
  if (!is.list(scdet_options)) {
    stop("'scdet_options' must be a list.", call. = FALSE)
  }
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
  if (!identical(max_mito, "SCdetMito")) {
    max_mito <- qc_bounds$max_mito
  }

  cutoff_plan <- build_multisample_cutoff_plan(
    seurat_obj = seurat_obj,
    sample_by = by,
    group_by = group_by,
    mitoRatio = mitoRatio,
      max_mito = max_mito,
      cutoff_strategy = cutoff_strategy,
    cutoff_quantile = cutoff_quantile,
    cutoff_support_fraction = cutoff_support_fraction,
    scdet_options = utils::modifyList(
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
    ),
    use_recommended_cutoff = use_recommended_cutoff,
    table_out = table_out,
    plot = plot && identical(max_mito, "SCdetMito"),
    output_dir = output_dir
  )
  applied_cutoffs <- cutoff_plan$sample_plan$applied_cutoff

  if (table_out) {
    utils::write.csv(
      cutoff_plan$sample_plan,
      file = file.path(output_dir, "SCQCmulti_cutoff_plan.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      cutoff_plan$group_plan,
      file = file.path(output_dir, "SCQCmulti_group_cutoffs.csv"),
      row.names = FALSE
    )
  }

  if (mode == "all") {
    reference_cell_count <- ncol(seurat_obj)
    cell_level_cutoffs <- build_cell_level_cutoff_vector(
      seurat_obj = seurat_obj,
      sample_by = by,
      sample_plan = cutoff_plan$sample_plan
    )
    before_summary <- build_qc_summary(
      seurat_obj,
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      mitoRatio = mitoRatio
    )
    before_summary$Stage <- "before"
    before_summary$Mode <- mode
    before_summary$MitoCutoff <- if (is.na(cutoff_plan$final_cutoff)) NA_real_ else cutoff_plan$final_cutoff
    before_summary <- add_qc_outcome_metrics(before_summary, reference_cell_count)
    before_summary <- add_cutoff_metrics(before_summary, applied_cutoffs, cutoff_strategy)

    if (plot) {
      plot_obj <- seurat_obj
      plot_obj@meta.data$.scdetmito_temp_group <- plot_obj@meta.data[[by]]
      SCQC_processedPlots(
        plot_obj,
        by = ".scdetmito_temp_group",
        flag = "M-Checked",
        nFeature_RNA = nFeature_RNA,
        nCount_RNA = nCount_RNA,
        mitoRatio = mitoRatio,
        output_dir = output_dir
      )
    }

    message("QC Summary [before]:")
    print(before_summary)
    message(
      "Applied mitoRatio cutoff range: ",
      paste(signif(range(applied_cutoffs), 3), collapse = " - ")
    )

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
      max_mito = max(applied_cutoffs, na.rm = TRUE),
      max_mito_by_cell = cell_level_cutoffs
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
    } else {
      doublet_status <- list(
        requested = FALSE,
        status = "not_requested",
        error_message = NA_character_
      )
    }

    after_summary <- build_qc_summary(
      filtered_obj,
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      mitoRatio = mitoRatio
    )
    after_summary$Stage <- "after"
    after_summary$Mode <- mode
    after_summary$MitoCutoff <- if (is.na(cutoff_plan$final_cutoff)) NA_real_ else cutoff_plan$final_cutoff
    after_summary <- add_qc_outcome_metrics(after_summary, reference_cell_count)
    after_summary <- add_cutoff_metrics(after_summary, applied_cutoffs, cutoff_strategy)

    if (plot) {
      plot_obj <- filtered_obj
      plot_obj@meta.data$.scdetmito_temp_group <- plot_obj@meta.data[[by]]
      SCQC_processedPlots(
        plot_obj,
        by = ".scdetmito_temp_group",
        flag = "M-Filtered",
        nFeature_RNA = nFeature_RNA,
        nCount_RNA = nCount_RNA,
        mitoRatio = mitoRatio,
        output_dir = output_dir
      )
    }

    message("QC Summary [after]:")
    print(after_summary)

    if (table_out) {
      utils::write.csv(
        rbind(before_summary, after_summary),
        file = file.path(output_dir, "SCQCmulti_all_summary.csv"),
        row.names = FALSE
      )
    }

    filtered_obj@misc$SCdetMito_QC <- attach_qc_provenance_to_cutoff_plan(
      cutoff_plan = cutoff_plan,
      provenance = build_qc_provenance(
        function_name = "SCQCmulti",
        input_cell_count = reference_cell_count,
        output_cell_count = ncol(filtered_obj),
        parameters = list(
          sample_col = by,
          group_col = group_by,
          mode = mode,
          cutoff_strategy = cutoff_strategy,
          cutoff_support_fraction = cutoff_support_fraction,
          feature_col = nFeature_RNA,
          count_col = nCount_RNA,
          mito_col = mitoRatio,
          min_genes = min_genes,
          max_genes = max_genes,
          min_counts = min_counts,
          max_counts = max_counts,
          min_mito = min_mito,
          max_mito = max_mito,
          use_recommended_cutoff = use_recommended_cutoff,
          remove_doublets = removeDouble
        ),
        cutoff_plan = cutoff_plan,
        cutoff_applied = cutoff_plan$cutoff_applied,
        cutoff_applied_source = cutoff_plan$cutoff_applied_source,
        cutoff_applied_strategy = cutoff_plan$cutoff_applied_strategy,
        applied_cutoff_level = cutoff_plan$applied_cutoff_level,
        final_selected_cutoff = cutoff_plan$final_selected_cutoff,
        final_recommended_cutoff = cutoff_plan$final_recommended_cutoff,
        cutoff_source = paste(unique(cutoff_plan$sample_plan$cutoff_source), collapse = ";"),
        cutoff_confidence = paste(unique(cutoff_plan$sample_plan$cutoff_confidence), collapse = ";"),
        loss_test = if (!is.null(cutoff_plan$detection)) cutoff_plan$detection$settings$loss_test %||% NA_character_ else NA_character_,
        sample_cutoff_method = if (!is.null(cutoff_plan$detection)) cutoff_plan$detection$settings$sample_cutoff_method %||% NA_character_ else NA_character_,
        fallback_events = cutoff_plan$sample_plan[cutoff_plan$sample_plan$fallback_used %in% TRUE, , drop = FALSE],
        doublet_filtering = doublet_status
      )
    )

    return(filtered_obj)
  }

  metadata <- seurat_obj@meta.data
  sample_levels <- unique(as.character(metadata[[by]]))

  if (plot) {
    plot_obj <- seurat_obj
    plot_obj@meta.data$.scdetmito_temp_group <- plot_obj@meta.data[[by]]
    SCQC_processedPlots(
      plot_obj,
      by = ".scdetmito_temp_group",
      flag = "M-Checked",
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      mitoRatio = mitoRatio,
      output_dir = output_dir
    )
  }

  split_results <- lapply(sample_levels, function(sample_name) {
    sample_cells <- rownames(metadata)[metadata[[by]] == sample_name]
    sample_obj <- subset(seurat_obj, cells = sample_cells)
    sample_cutoff <- cutoff_plan$sample_plan$applied_cutoff[
      match(sample_name, cutoff_plan$sample_plan$sample_id)
    ]
    SCQCone(
      sample_obj,
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      mitoRatio = mitoRatio,
      min_genes = min_genes,
      max_genes = max_genes,
      min_counts = min_counts,
      max_counts = max_counts,
      min_mito = min_mito,
      max_mito = sample_cutoff,
      removeDouble = removeDouble,
      plot = FALSE,
      table_out = FALSE,
      output_dir = output_dir,
      fail_action = fail_action
    )
  })
  names(split_results) <- sample_levels

  qcpassed_seurat_obj <- if (length(split_results) == 1) {
    split_results[[1]]
  } else {
    Reduce(function(x, y) merge(x = x, y = y), split_results)
  }

  if (plot) {
    plot_obj <- qcpassed_seurat_obj
    plot_obj@meta.data$.scdetmito_temp_group <- plot_obj@meta.data[[by]]
    SCQC_processedPlots(
      plot_obj,
      by = ".scdetmito_temp_group",
      flag = "M-Filtered",
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      mitoRatio = mitoRatio,
      output_dir = output_dir
    )
  }

  if (table_out) {
    summary_table <- data.frame(
      sample = sample_levels,
      group = cutoff_plan$sample_plan$group_id[match(sample_levels, cutoff_plan$sample_plan$sample_id)],
      applied_cutoff = cutoff_plan$sample_plan$applied_cutoff[match(sample_levels, cutoff_plan$sample_plan$sample_id)],
      cutoff_strategy = cutoff_strategy,
      before_cells = vapply(sample_levels, function(sample_name) {
        sum(metadata[[by]] == sample_name)
      }, numeric(1)),
      after_cells = vapply(sample_levels, function(sample_name) {
        sum(qcpassed_seurat_obj@meta.data[[by]] == sample_name, na.rm = TRUE)
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
    summary_table$cells_removed <- summary_table$before_cells - summary_table$after_cells
    summary_table$retention_rate <- round(
      summary_table$after_cells / summary_table$before_cells,
      4
    )
    utils::write.csv(
      summary_table,
      file = file.path(output_dir, "SCQCmulti_split_summary.csv"),
      row.names = FALSE
    )
  }

  split_doublet_status <- lapply(split_results, function(obj) {
    obj@misc$SCdetMito_QC$doublet_filtering %||% NULL
  })
  qcpassed_seurat_obj@misc$SCdetMito_QC <- attach_qc_provenance_to_cutoff_plan(
    cutoff_plan = cutoff_plan,
    provenance = build_qc_provenance(
      function_name = "SCQCmulti",
      input_cell_count = ncol(seurat_obj),
      output_cell_count = ncol(qcpassed_seurat_obj),
      parameters = list(
        sample_col = by,
        group_col = group_by,
        mode = mode,
        cutoff_strategy = cutoff_strategy,
        cutoff_support_fraction = cutoff_support_fraction,
        feature_col = nFeature_RNA,
        count_col = nCount_RNA,
        mito_col = mitoRatio,
        min_genes = min_genes,
        max_genes = max_genes,
        min_counts = min_counts,
        max_counts = max_counts,
        min_mito = min_mito,
        max_mito = max_mito,
        use_recommended_cutoff = use_recommended_cutoff,
        remove_doublets = removeDouble
      ),
      cutoff_plan = cutoff_plan,
      cutoff_applied = cutoff_plan$cutoff_applied,
      cutoff_applied_source = cutoff_plan$cutoff_applied_source,
      cutoff_applied_strategy = cutoff_plan$cutoff_applied_strategy,
      applied_cutoff_level = cutoff_plan$applied_cutoff_level,
      final_selected_cutoff = cutoff_plan$final_selected_cutoff,
      final_recommended_cutoff = cutoff_plan$final_recommended_cutoff,
      cutoff_source = paste(unique(cutoff_plan$sample_plan$cutoff_source), collapse = ";"),
      cutoff_confidence = paste(unique(cutoff_plan$sample_plan$cutoff_confidence), collapse = ";"),
      loss_test = if (!is.null(cutoff_plan$detection)) cutoff_plan$detection$settings$loss_test %||% NA_character_ else NA_character_,
      sample_cutoff_method = if (!is.null(cutoff_plan$detection)) cutoff_plan$detection$settings$sample_cutoff_method %||% NA_character_ else NA_character_,
      fallback_events = cutoff_plan$sample_plan[cutoff_plan$sample_plan$fallback_used %in% TRUE, , drop = FALSE],
      doublet_filtering = list(
        requested = isTRUE(removeDouble),
        status = if (isTRUE(removeDouble)) "per_sample" else "not_requested",
        per_sample = split_doublet_status
      )
    )
  )

  qcpassed_seurat_obj
}
