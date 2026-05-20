#' SCQC_processedPlots: draw QC plots for processed Seurat objects
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param by Metadata column used to group cells in the plots.
#' @param flag Prefix used for output file names.
#' @param nFeature_RNA Metadata column storing detected feature counts.
#' @param nCount_RNA Metadata column storing UMI/count totals.
#' @param mitoRatio Metadata column storing mitochondrial ratios.
#' @param output_dir Directory for exported plots. Defaults to `"."`.
#'
#' @return Invisibly returns the generated `ggplot2` plot objects.
#' @export
#'
#' @examples
#' # DO NOT RUN
#' # library(Seurat)
#' # seu <- CreateSeuratObject(matrix(rpois(2000, 5), nrow = 100))
#' # seu$sample <- rep(c("A", "B"), each = ncol(seu) / 2)
#' # seu$mitoRatio <- runif(ncol(seu))
#' # SCQC_processedPlots(seu, by = "sample", flag = "demo", output_dir = tempdir())
#' # DO NOT RUN
SCQC_processedPlots <- function(seurat_obj,
                                by,
                                flag,
                                nFeature_RNA = "nFeature_RNA",
                                nCount_RNA = "nCount_RNA",
                                mitoRatio = "mitoRatio",
                                output_dir = ".") {
  output_dir <- ensure_output_dir(output_dir)
  seurat_obj <- check_seu(seurat_obj, by)
  seurat_obj <- check_seu(seurat_obj, nFeature_RNA, must_be_numeric = TRUE)
  seurat_obj <- check_seu(seurat_obj, nCount_RNA, must_be_numeric = TRUE)
  seurat_obj <- check_seu(
    seurat_obj,
    mitoRatio,
    normalize_fraction = TRUE,
    must_be_numeric = TRUE
  )

  metadata <- data.frame(
    sample_group = as.character(meta_column(seurat_obj, by)),
    plot_nCount_RNA = meta_column(seurat_obj, nCount_RNA),
    plot_nFeature_RNA = meta_column(seurat_obj, nFeature_RNA),
    plot_mitoRatio = meta_column(seurat_obj, mitoRatio),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  metadata$log10GenesPerUMI <- NA_real_
  valid_complexity <- metadata$plot_nFeature_RNA > 0 & metadata$plot_nCount_RNA > 1
  metadata$log10GenesPerUMI[valid_complexity] <-
    log10(metadata$plot_nFeature_RNA[valid_complexity]) /
      log10(metadata$plot_nCount_RNA[valid_complexity])

  sample_levels <- unique(metadata$sample_group)
  cols <- scdetmito_palette(length(sample_levels))
  names(cols) <- sample_levels

  p_cells <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(x = .data$sample_group, fill = .data$sample_group)
  ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::geom_bar(width = 0.6) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::ggtitle("Cells number of each sample") +
    Seurat::NoLegend() +
    ggplot2::labs(x = "Samples", y = "Counts")

  p_umi_density <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(x = .data$plot_nCount_RNA, fill = .data$sample_group)
  ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::geom_density(alpha = 0.67) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
    ggplot2::scale_x_log10() +
    ggplot2::ggtitle("Cell density per UMI number") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Samples")) +
    ggplot2::labs(x = "UMI numbers", y = "Cell density")

  p_gene_density <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(x = .data$plot_nFeature_RNA, fill = .data$sample_group)
  ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::geom_density(alpha = 0.67) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
    ggplot2::scale_x_log10() +
    ggplot2::ggtitle("Cell density per gene number") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Samples")) +
    ggplot2::labs(x = "Gene numbers", y = "Cell density")

  p_genes_per_cell <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(
      x = .data$sample_group,
      y = log10(.data$plot_nFeature_RNA),
      fill = .data$sample_group
    )
  ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::ggtitle("Genes per cell in each sample") +
    ggplot2::labs(x = "Samples", y = "log10(Gene number)")

  p_gene_vs_umi <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(
      x = .data$plot_nCount_RNA,
      y = .data$plot_nFeature_RNA,
      color = .data$plot_mitoRatio
    )
  ) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::scale_color_gradient(low = "lightgray", high = "tomato") +
    ggplot2::stat_smooth(method = stats::lm, se = FALSE) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
    ggplot2::facet_wrap(~sample_group)

  p_mito <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(x = .data$plot_mitoRatio, fill = .data$sample_group)
  ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::geom_density(alpha = 0.5, linewidth = 0.1) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Samples")) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::ggtitle("Cells distribution by MitoRatio") +
    ggplot2::labs(x = "MitoRatio", y = "Density")

  p_complexity <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(x = .data$log10GenesPerUMI, fill = .data$sample_group)
  ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
    ggplot2::geom_density(alpha = 0.5, linewidth = 0.2) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::ggtitle("Genes complexity") +
    ggplot2::labs(x = "log10(Genes Per UMI)", y = "Density") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Samples"))

  plots <- list(
    cells = p_cells,
    umi_density = p_umi_density,
    gene_density = p_gene_density,
    genes_per_cell = p_genes_per_cell,
    gene_vs_umi = p_gene_vs_umi,
    mito = p_mito,
    complexity = p_complexity
  )

  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(flag, "_Ncells_in_each_sample.pdf")),
    plot = p_cells,
    width = 8,
    height = 5
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(flag, "_Cell_density_per_UMI.pdf")),
    plot = p_umi_density,
    width = 8,
    height = 5
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(flag, "_Cell_density_per_Gene.pdf")),
    plot = p_gene_density,
    width = 8,
    height = 5
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(flag, "_gene_per_cell.pdf")),
    plot = p_genes_per_cell,
    width = 8,
    height = 5
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(flag, "_gene_vs_umi-mitoRatio.pdf")),
    plot = p_gene_vs_umi,
    width = 12,
    height = 10
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(flag, "_mitochondrial.pdf")),
    plot = p_mito,
    width = 9,
    height = 5
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(flag, "_genes_Complexity.pdf")),
    plot = p_complexity,
    width = 9,
    height = 5
  )

  invisible(plots)
}
