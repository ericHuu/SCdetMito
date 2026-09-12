# Maintainer script for building the bundled SCdetMito PBMC demo.
#
# The output is a downsampled multi-donor Seurat object derived from public
# 10x Genomics GEM-X Human PBMC donor datasets. This script is not run during
# package installation, vignettes, tests, or R CMD check.

build_demo_pbmc_gemx_multidonor <- function(
  output_file = file.path("inst", "extdata", "demo_pbmc_multidonor_seurat.rds"),
  cache_dir = file.path("data-public", "gemx_pbmc"),
  local_files = NULL,
  downsample_per_donor = 250,
  seed = 1401,
  min_cells = 0,
  min_features = 0,
  verbose = TRUE
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required to build the demo PBMC object.", call. = FALSE)
  }

  if (!exists("ensure_mito_ratio", mode = "function")) {
    if (requireNamespace("devtools", quietly = TRUE)) {
      devtools::load_all(".", quiet = TRUE)
    } else {
      stop("Run this script from the package root or install SCdetMito first.", call. = FALSE)
    }
  }

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

  donor_info <- data.frame(
    donor = paste0("donor", 1:4),
    sample = paste0("donor", 1:4),
    sex = c("male", "male", "female", "female"),
    age_group = c("18_35", "18_35", "36_50", "36_50"),
    stringsAsFactors = FALSE
  )
  donor_info$group <- paste(donor_info$sex, donor_info$age_group, sep = "_")

  dataset_ids <- paste0(
    "5k_Human_Donor",
    1:4,
    "_PBMC_3p_gem-x_5k_Human_Donor",
    1:4,
    "_PBMC_3p_gem-x"
  )
  default_urls <- paste0(
    "https://cf.10xgenomics.com/samples/cell-exp/9.0.0/",
    dataset_ids,
    "/",
    dataset_ids,
    "_count_sample_filtered_feature_bc_matrix.h5"
  )
  default_local_files <- file.path(
    cache_dir,
    paste0(donor_info$donor, "_filtered_feature_bc_matrix.h5")
  )

  if (is.null(local_files)) {
    local_files <- default_local_files
  }
  if (length(local_files) != nrow(donor_info)) {
    stop("'local_files' must contain one path per donor.", call. = FALSE)
  }

  seurat_list <- vector("list", nrow(donor_info))
  for (i in seq_len(nrow(donor_info))) {
    donor <- donor_info$donor[i]
    local_file <- local_files[i]
    if (!file.exists(local_file)) {
      if (isTRUE(verbose)) {
        message("Downloading ", donor, " filtered feature-barcode matrix.")
      }
      status <- tryCatch(
        utils::download.file(default_urls[i], local_file, mode = "wb", quiet = !verbose),
        error = function(e) e
      )
      if (inherits(status, "error") || !file.exists(local_file)) {
        stop(
          "Could not download the public 10x matrix for ",
          donor,
          ". URL attempted: ",
          default_urls[i],
          ". Provide the file manually through 'local_files'.",
          call. = FALSE
        )
      }
    } else if (isTRUE(verbose)) {
      message("Using cached file for ", donor, ": ", local_file)
    }

    counts <- Seurat::Read10X_h5(local_file)
    if (is.list(counts)) {
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else {
        counts <- counts[[1]]
      }
    }

    obj <- Seurat::CreateSeuratObject(
      counts = counts,
      project = donor,
      min.cells = min_cells,
      min.features = min_features
    )
    obj$sample <- donor_info$sample[i]
    obj$donor <- donor_info$donor[i]
    obj$sex <- donor_info$sex[i]
    obj$age_group <- donor_info$age_group[i]
    obj$group <- donor_info$group[i]
    obj$condition <- "healthy_pbmc"
    obj$tissue <- "PBMC"
    obj$species <- "human"
    obj$dataset <- "10x_gemx_pbmc_donor1_4_demo"
    obj$data_source <- "10x Genomics public GEM-X PBMC donor datasets"
    obj$demo_note <- "for software demonstration and testing only"

    if (!is.null(downsample_per_donor) && ncol(obj) > downsample_per_donor) {
      set.seed(seed + i)
      keep_cells <- sample(colnames(obj), downsample_per_donor)
      obj <- subset(obj, cells = keep_cells)
    }
    seurat_list[[i]] <- obj
  }

  merged <- Reduce(
    function(x, y) merge(x, y = y, merge.data = FALSE),
    seurat_list
  )
  merged <- ensure_mito_ratio(
    merged,
    mito_col = "mitoRatio",
    species = "human",
    overwrite = TRUE,
    verbose = verbose
  )
  merged@misc$SCdetMito_demo_provenance <- list(
    source = "10x Genomics public GEM-X Human PBMC Donor 1-4 datasets",
    license = "CC BY 4.0",
    cell_ranger_version = "9.0.0",
    downsample_per_donor = downsample_per_donor,
    seed = seed,
    created = as.character(Sys.time()),
    urls = stats::setNames(default_urls, donor_info$donor)
  )

  saveRDS(merged, output_file, compress = "xz")
  if (isTRUE(verbose)) {
    message("Saved demo object: ", output_file)
    message("Cells: ", ncol(merged), "; features: ", nrow(merged))
  }
  invisible(merged)
}

if (identical(environment(), globalenv()) && !interactive()) {
  build_demo_pbmc_gemx_multidonor()
}
