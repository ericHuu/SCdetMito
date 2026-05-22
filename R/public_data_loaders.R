# SCdetMito
# Author: Silu Hu
# Contact: husilu0902@gmail.com
# Version: 1.4.2
# Last updated: 2026-05-22

# Optional public data loaders for SCdetMito examples. Online datasets are
# downloaded only when users explicitly call these helpers.

#' Load the public PBMC3K dataset for SCdetMito examples
#'
#' @description
#' This optional helper loads the public 10x Genomics PBMC3K dataset either
#' through SeuratData or by downloading the filtered gene-barcode matrix from
#' the public 10x Genomics URL. PBMC3K is provided as a small public
#' smoke-test example. Because PBMCs often have relatively low mitochondrial
#' fractions and dataset-specific QC structure, PBMC3K cutoff results should be
#' interpreted as a functionality check rather than evidence of a universal
#' mitochondrial threshold.
#'
#' @param method Loading method. `"SeuratData"` loads an already installed
#'   SeuratData PBMC3K object. `"10x"` downloads the public 10x Genomics
#'   filtered gene-barcode matrix to `cache_dir`.
#' @param cache_dir Directory used to cache the downloaded 10x archive and
#'   extracted matrix. Defaults to the user cache directory for SCdetMito.
#' @param add_mito Whether to add mitochondrial ratio metadata when `mito_col`
#'   is absent. Defaults to `TRUE`.
#' @param mito_col Metadata column used for mitochondrial ratios. Defaults to
#'   `"mitoRatio"`.
#' @param min_cells Minimum cells per feature passed to
#'   [Seurat::CreateSeuratObject()] for the 10x method.
#' @param min_features Minimum detected features per cell passed to
#'   [Seurat::CreateSeuratObject()] for the 10x method.
#' @param verbose Whether to print cache and download messages. Defaults to
#'   `TRUE`.
#'
#' @return A Seurat object containing the public PBMC3K data with minimal
#'   metadata and mitochondrial metadata when `add_mito = TRUE`.
#'
#' @details
#' This function is intended for explicit real-data smoke tests. It is not used
#' during package installation, examples, automated tests, or R CMD check.
#' Internet access is required for `method = "10x"` unless the 10x archive has
#' already been cached.
#'
#' @examples
#' \dontrun{
#' pbmc <- load_pbmc3k_online(method = "10x")
#' det <- SCdetMito(pbmc, sample_col = "sample", mito_col = "mitoRatio")
#' }
#' @export
load_pbmc3k_online <- function(method = c("SeuratData", "10x"),
                               cache_dir = tools::R_user_dir("SCdetMito", which = "cache"),
                               add_mito = TRUE,
                               mito_col = "mitoRatio",
                               min_cells = 3,
                               min_features = 200,
                               verbose = TRUE) {
  method <- match.arg(method)
  .ensure_seurat_available()

  seurat_obj <- switch(method,
    SeuratData = .load_seuratdata_dataset(
      dataset = "pbmc3k",
      missing_package_message = paste0(
        "SeuratData is required for method = 'SeuratData'. Install it with ",
        "remotes::install_github('satijalab/seurat-data'), then run ",
        "SeuratData::InstallData('pbmc3k')."
      ),
      missing_dataset_message = "The pbmc3k dataset is not installed. Run SeuratData::InstallData('pbmc3k') and retry."
    ),
    `10x` = .download_and_read_10x_matrix(
      urls = "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
      dataset_id = "pbmc3k",
      cache_dir = cache_dir,
      project = "pbmc3k",
      min_cells = min_cells,
      min_features = min_features,
      verbose = verbose
    )
  )

  seurat_obj$sample <- "pbmc3k"
  seurat_obj$condition <- "public_pbmc"
  seurat_obj$dataset <- "10x_pbmc3k"
  seurat_obj$data_source <- if (identical(method, "SeuratData")) {
    "SeuratData pbmc3k"
  } else {
    "10x Genomics public PBMC3K"
  }

  .add_or_normalize_mito_ratio(
    seurat_obj = seurat_obj,
    add_mito = add_mito,
    mito_col = mito_col,
    species = "human",
    dataset_label = "PBMC3K"
  )
}

#' Load a public 10x Genomics mouse heart 10k dataset
#'
#' @description
#' This optional helper loads the public 10x Genomics "10k Heart Cells from an
#' E18 mouse" dataset. The dataset is not included in SCdetMito and is
#' downloaded only when this function is called. It is intended as an optional
#' solid-tissue stress-test example because this dataset may have high
#' mitochondrial burden and stringent QC behavior.
#'
#' @param cache_dir Directory used to cache the downloaded 10x archive and
#'   extracted matrix. Defaults to the user cache directory for SCdetMito.
#' @param add_mito Whether to add mitochondrial ratio metadata when `mito_col`
#'   is absent. Defaults to `TRUE`.
#' @param mito_col Metadata column used for mitochondrial ratios. Defaults to
#'   `"mitoRatio"`.
#' @param min_cells Minimum cells per feature passed to
#'   [Seurat::CreateSeuratObject()].
#' @param min_features Minimum detected features per cell passed to
#'   [Seurat::CreateSeuratObject()].
#' @param url Optional URL or vector of URLs overriding the built-in public
#'   10x Genomics URL.
#' @param local_file Optional path to a manually downloaded filtered
#'   feature-barcode matrix archive. When supplied, no download is attempted.
#' @param verbose Whether to print cache and download messages. Defaults to
#'   `TRUE`.
#'
#' @return A Seurat object containing public mouse heart cells with
#'   mitochondrial metadata when `add_mito = TRUE`.
#'
#' @details
#' Internet access is required unless the 10x archive has already been cached.
#' The downloaded archive and extracted matrix are stored under `cache_dir`,
#' not inside the installed SCdetMito package.
#'
#' @examples
#' \dontrun{
#' heart <- load_heart10k_online()
#' det <- SCdetMito(heart, sample_col = "sample", mito_col = "mitoRatio")
#' }
#' @export
load_heart10k_online <- function(cache_dir = tools::R_user_dir("SCdetMito", which = "cache"),
                                 add_mito = TRUE,
                                 mito_col = "mitoRatio",
                                 min_cells = 3,
                                 min_features = 200,
                                 url = NULL,
                                 local_file = NULL,
                                 verbose = TRUE) {
  .ensure_seurat_available()

  seurat_obj <- .download_and_read_10x_matrix(
    urls = url %||% c(
      "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_10k_v3/heart_10k_v3_filtered_feature_bc_matrix.tar.gz",
      "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_10k_v3/heart_10k_v3_filtered_feature_bc_matrix.tar.gz"
    ),
    dataset_id = "heart10k",
    cache_dir = cache_dir,
    project = "heart10k",
    min_cells = min_cells,
    min_features = min_features,
    local_file = local_file,
    verbose = verbose
  )

  seurat_obj$sample <- "heart10k"
  seurat_obj$condition <- "public_mouse_heart"
  seurat_obj$dataset <- "10x_heart10k_e18_mouse"
  seurat_obj$tissue <- "heart"
  seurat_obj$species <- "mouse"
  seurat_obj$data_source <- "10x Genomics public Heart10k"

  .add_or_normalize_mito_ratio(
    seurat_obj = seurat_obj,
    add_mito = add_mito,
    mito_col = mito_col,
    species = "mouse",
    dataset_label = "10x Genomics mouse heart 10k"
  )
}

#' Load a public 10x Genomics mouse brain 5k nuclei dataset
#'
#' @description
#' This optional helper loads the public 10x Genomics 5k adult mouse brain
#' nuclei dataset isolated with the Chromium Nuclei Isolation Kit. The dataset
#' is not included in SCdetMito and is downloaded only when this function is
#' called. It is intended as a real tissue/nuclei candidate for dataset
#' scouting before selecting SoftwareX figures.
#'
#' @param cache_dir Directory used to cache the downloaded 10x HDF5 matrix.
#'   Defaults to the user cache directory for SCdetMito.
#' @param add_mito Whether to add mitochondrial ratio metadata when `mito_col`
#'   is absent. Defaults to `TRUE`.
#' @param mito_col Metadata column used for mitochondrial ratios. Defaults to
#'   `"mitoRatio"`.
#' @param min_cells Minimum cells per feature passed to
#'   [Seurat::CreateSeuratObject()].
#' @param min_features Minimum detected features per cell passed to
#'   [Seurat::CreateSeuratObject()].
#' @param url Optional URL overriding the built-in public 10x Genomics HDF5 URL.
#' @param local_file Optional path to a manually downloaded filtered
#'   feature-barcode HDF5 file. When supplied, no download is attempted.
#' @param verbose Whether to print cache and download messages. Defaults to
#'   `TRUE`.
#'
#' @return A Seurat object containing public mouse brain nuclei with
#'   mitochondrial metadata when `add_mito = TRUE`.
#'
#' @details
#' Internet access is required unless the 10x HDF5 file has already been
#' cached. Nuclei datasets may have low mitochondrial fractions, so this loader
#' should be treated as a dataset-scouting candidate rather than a guaranteed
#' high-mitochondrial example.
#'
#' @examples
#' \dontrun{
#' brain <- load_mousebrain5k_online()
#' det <- SCdetMito(brain, sample_col = "sample", mito_col = "mitoRatio")
#' }
#' @export
load_mousebrain5k_online <- function(cache_dir = tools::R_user_dir("SCdetMito", which = "cache"),
                                     add_mito = TRUE,
                                     mito_col = "mitoRatio",
                                     min_cells = 3,
                                     min_features = 200,
                                     url = NULL,
                                     local_file = NULL,
                                     verbose = TRUE) {
  .ensure_seurat_available()

  seurat_obj <- .download_and_read_10x_h5(
    urls = url %||% "https://cf.10xgenomics.com/samples/cell-exp/7.0.0/5k_mouse_brain_CNIK_3pv3/5k_mouse_brain_CNIK_3pv3_filtered_feature_bc_matrix.h5",
    dataset_id = "mousebrain5k",
    cache_dir = cache_dir,
    project = "mousebrain5k",
    min_cells = min_cells,
    min_features = min_features,
    local_file = local_file,
    verbose = verbose
  )

  seurat_obj$sample <- "mousebrain5k"
  seurat_obj$condition <- "public_mouse_brain"
  seurat_obj$dataset <- "10x_mousebrain5k"
  seurat_obj$tissue <- "brain"
  seurat_obj$species <- "mouse"
  seurat_obj$data_source <- "10x Genomics public mouse brain 5k"

  .add_or_normalize_mito_ratio(
    seurat_obj = seurat_obj,
    add_mito = add_mito,
    mito_col = mito_col,
    species = "mouse",
    dataset_label = "10x Genomics mouse brain 5k"
  )
}

.ensure_seurat_available <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required to load public data as a Seurat object.", call. = FALSE)
  }
}

.load_seuratdata_dataset <- function(dataset,
                                     missing_package_message,
                                     missing_dataset_message,
                                     install_if_missing = FALSE,
                                     verbose = TRUE) {
  if (!requireNamespace("SeuratData", quietly = TRUE)) {
    stop(missing_package_message, call. = FALSE)
  }

  if (!.seuratdata_dataset_installed(dataset)) {
    if (isTRUE(install_if_missing)) {
      if (isTRUE(verbose)) {
        message("Installing SeuratData dataset: ", dataset)
      }
      install_status <- tryCatch(
        getExportedValue("SeuratData", "InstallData")(dataset),
        error = function(e) e
      )
      if (inherits(install_status, "error") || !.seuratdata_dataset_installed(dataset)) {
        stop(
          "Failed to install SeuratData dataset '",
          dataset,
          "'. ",
          if (inherits(install_status, "error")) conditionMessage(install_status) else missing_dataset_message,
          call. = FALSE
        )
      }
    } else {
      stop(missing_dataset_message, call. = FALSE)
    }
  }

  tryCatch(
    getExportedValue("SeuratData", "LoadData")(dataset),
    error = function(e) {
      stop(
        "Failed to load SeuratData dataset '",
        dataset,
        "': ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

.seuratdata_dataset_installed <- function(dataset) {
  installed_data <- tryCatch(
    getExportedValue("SeuratData", "InstalledData")(),
    error = function(e) NULL
  )
  installed_values <- unique(as.character(unlist(installed_data, use.names = FALSE)))
  any(grepl(paste0("^", dataset), installed_values, ignore.case = TRUE))
}

.download_and_read_10x_matrix <- function(urls,
                                          dataset_id,
                                          cache_dir,
                                          project,
                                          min_cells,
                                          min_features,
                                          local_file = NULL,
                                          verbose) {
  .ensure_seurat_available()
  .validate_cache_dir(cache_dir)

  extract_dir <- file.path(cache_dir, dataset_id)
  if (!dir.exists(extract_dir)) {
    dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(extract_dir)) {
    stop("Could not create extraction directory for ", dataset_id, ": ", extract_dir, call. = FALSE)
  }

  archive_file <- if (!is.null(local_file)) {
    .validate_local_public_data_file(local_file, dataset_id)
  } else {
    file.path(cache_dir, basename(urls[[1]]))
  }
  if (is.null(local_file) && (!file.exists(archive_file) || isTRUE(file.info(archive_file)$size <= 0))) {
    .download_archive_with_fallback(
      urls = urls,
      destfile = archive_file,
      dataset_id = dataset_id,
      verbose = verbose
    )
  } else if (!is.null(local_file) && isTRUE(verbose)) {
    message("Using local file for ", dataset_id, ": ", archive_file)
  } else if (isTRUE(verbose)) {
    message("Using cached file for ", dataset_id, ": ", archive_file)
  }

  matrix_dir <- .find_10x_matrix_dir(extract_dir)
  if (is.null(matrix_dir)) {
    if (isTRUE(verbose)) {
      message("Extracting ", dataset_id, " archive under: ", extract_dir)
    }
    extract_status <- tryCatch(
      utils::untar(archive_file, exdir = extract_dir),
      error = function(e) e
    )
    if (inherits(extract_status, "error")) {
      stop(
        "Failed to extract ",
        dataset_id,
        " archive: ",
        conditionMessage(extract_status),
        call. = FALSE
      )
    }
    matrix_dir <- .find_10x_matrix_dir(extract_dir)
  }
  if (is.null(matrix_dir)) {
    stop(
      "Could not locate the extracted 10x matrix directory for ",
      dataset_id,
      " under: ",
      extract_dir,
      call. = FALSE
    )
  }
  if (isTRUE(verbose)) {
    message("Matrix directory detected for ", dataset_id, ": ", matrix_dir)
  }

  counts <- tryCatch(
    Seurat::Read10X(data.dir = matrix_dir),
    error = function(e) {
      stop("Failed to read 10x matrix for ", dataset_id, ": ", conditionMessage(e), call. = FALSE)
    }
  )
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  }

  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = project,
    min.cells = min_cells,
    min.features = min_features
  )
  if (isTRUE(verbose)) {
    message(
      "Loaded ",
      ncol(seurat_obj),
      " cells and ",
      nrow(seurat_obj),
      " features for ",
      dataset_id,
      "."
    )
  }
  seurat_obj
}

.download_and_read_10x_h5 <- function(urls,
                                      dataset_id,
                                      cache_dir,
                                      project,
                                      min_cells,
                                      min_features,
                                      local_file = NULL,
                                      verbose) {
  .ensure_seurat_available()
  .validate_cache_dir(cache_dir)

  h5_file <- if (!is.null(local_file)) {
    .validate_local_public_data_file(local_file, dataset_id)
  } else {
    file.path(cache_dir, basename(urls[[1]]))
  }
  if (is.null(local_file) && (!file.exists(h5_file) || isTRUE(file.info(h5_file)$size <= 0))) {
    .download_archive_with_fallback(
      urls = urls,
      destfile = h5_file,
      dataset_id = dataset_id,
      verbose = verbose
    )
  } else if (!is.null(local_file) && isTRUE(verbose)) {
    message("Using local file for ", dataset_id, ": ", h5_file)
  } else if (isTRUE(verbose)) {
    message("Using cached file for ", dataset_id, ": ", h5_file)
  }

  counts <- tryCatch(
    Seurat::Read10X_h5(filename = h5_file),
    error = function(e) {
      stop("Failed to read 10x HDF5 matrix for ", dataset_id, ": ", conditionMessage(e), call. = FALSE)
    }
  )
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  }

  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = project,
    min.cells = min_cells,
    min.features = min_features
  )
  if (isTRUE(verbose)) {
    message(
      "Loaded ",
      ncol(seurat_obj),
      " cells and ",
      nrow(seurat_obj),
      " features for ",
      dataset_id,
      "."
    )
  }
  seurat_obj
}

.download_archive_with_fallback <- function(urls,
                                            destfile,
                                            dataset_id,
                                            verbose) {
  attempted <- character(0)
  for (idx in seq_along(urls)) {
    current_url <- urls[[idx]]
    attempted <- c(attempted, current_url)
    if (isTRUE(verbose)) {
      if (idx == 1L) {
        message("Downloading ", dataset_id, " from: ", current_url)
      } else {
        message("Trying fallback URL for ", dataset_id, ": ", current_url)
      }
    }
    status <- tryCatch(
      suppressWarnings(utils::download.file(
        url = current_url,
        destfile = destfile,
        mode = "wb",
        quiet = !isTRUE(verbose)
      )),
      error = function(e) e
    )
    if (!inherits(status, "error") &&
      identical(as.integer(status), 0L) &&
      file.exists(destfile) &&
      isTRUE(file.info(destfile)$size > 0)) {
      return(invisible(destfile))
    }
    if (file.exists(destfile) && isTRUE(file.info(destfile)$size <= 0)) {
      unlink(destfile)
    }
  }

  stop(
    "Failed to download ",
    dataset_id,
    ". URLs attempted: ",
    paste(attempted, collapse = "; "),
    ". Check internet access or manually download the 10x filtered feature-barcode matrix and pass it with 'local_file'.",
    call. = FALSE
  )
}

.validate_local_public_data_file <- function(local_file, dataset_id) {
  if (!is.character(local_file) || length(local_file) != 1 || !nzchar(local_file)) {
    stop("'local_file' must be a single non-empty file path for ", dataset_id, ".", call. = FALSE)
  }
  if (!file.exists(local_file)) {
    stop("Local file for ", dataset_id, " does not exist: ", local_file, call. = FALSE)
  }
  if (!isTRUE(file.info(local_file)$size > 0)) {
    stop("Local file for ", dataset_id, " is empty: ", local_file, call. = FALSE)
  }
  normalizePath(local_file, winslash = "/", mustWork = TRUE)
}

.validate_cache_dir <- function(cache_dir) {
  if (!is.character(cache_dir) || length(cache_dir) != 1 || !nzchar(cache_dir)) {
    stop("'cache_dir' must be a non-empty directory path.", call. = FALSE)
  }
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(cache_dir)) {
    stop("Could not create cache directory: ", cache_dir, call. = FALSE)
  }
}

.find_10x_matrix_dir <- function(root_dir) {
  candidates <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
  candidates <- candidates[vapply(candidates, .is_10x_matrix_dir, logical(1))]
  if (!length(candidates)) {
    return(NULL)
  }
  candidates[[1]]
}

.is_10x_matrix_dir <- function(path) {
  has_matrix <- file.exists(file.path(path, "matrix.mtx")) ||
    file.exists(file.path(path, "matrix.mtx.gz"))
  has_barcodes <- file.exists(file.path(path, "barcodes.tsv")) ||
    file.exists(file.path(path, "barcodes.tsv.gz"))
  has_features <- file.exists(file.path(path, "genes.tsv")) ||
    file.exists(file.path(path, "genes.tsv.gz")) ||
    file.exists(file.path(path, "features.tsv")) ||
    file.exists(file.path(path, "features.tsv.gz"))

  has_matrix && has_barcodes && has_features
}

.add_or_normalize_mito_ratio <- function(seurat_obj,
                                         add_mito,
                                         mito_col,
                                         species,
                                         dataset_label) {
  if (!isTRUE(add_mito)) {
    return(seurat_obj)
  }

  prepared <- tryCatch(
    ensure_mito_ratio(
      object = seurat_obj,
      mito_col = mito_col,
      species = species,
      overwrite = FALSE,
      convert_percent = TRUE,
      verbose = FALSE
    ),
    error = function(e) e
  )
  if (!inherits(prepared, "error")) {
    return(prepared)
  }

  stop(
    "Failed to prepare mitochondrial ratio for ",
    dataset_label,
    ": ",
    conditionMessage(prepared),
    call. = FALSE
  )
}
