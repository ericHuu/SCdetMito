# SCdetMito
# Author: Silu Hu
# Contact: husilu0902@gmail.com
# Version: 1.4.3
# Last updated: 2026-05-23

common_mito_patterns <- function() {
  c("^MT-", "^mt-", "^Mt-", "^MT\\.", "^MT_", "^mt\\.", "^mt_", "^mitochondrial")
}

# Exact presets are fallback lists. Pattern-based detection remains the
# preferred default and may capture additional mitochondrial features retained
# in the object.
mito_species_exact_features <- function(species) {
  species <- match.arg(
    species %||% "auto",
    c(
      "auto", "human", "mouse", "rat", "pig", "cow", "cattle",
      "sheep", "goat", "chicken", "duck", "goose", "mixed"
    )
  )

  bare_mt_symbols <- c(
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6",
    "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB",
    "12S", "16S", "RNR1", "RNR2"
  )

  human_symbols <- c(
    "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8",
    "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4",
    "MT-ND5", "MT-ND6", "MT-CYB", "MT-RNR1", "MT-RNR2"
  )
  mouse_symbols <- c(
    "mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp8",
    "mt-Atp6", "mt-Co3", "mt-Nd3", "mt-Nd4l", "mt-Nd4",
    "mt-Nd5", "mt-Nd6", "mt-Cytb", "mt-Rnr1", "mt-Rnr2"
  )

  switch(species,
    human = human_symbols,
    mouse = unique(c(mouse_symbols, human_symbols)),
    rat = unique(c(mouse_symbols, human_symbols)),
    pig = bare_mt_symbols,
    cow = bare_mt_symbols,
    cattle = bare_mt_symbols,
    sheep = bare_mt_symbols,
    goat = bare_mt_symbols,
    chicken = bare_mt_symbols,
    duck = bare_mt_symbols,
    goose = bare_mt_symbols,
    mixed = unique(c(bare_mt_symbols, human_symbols, mouse_symbols)),
    auto = character(0)
  )
}

mito_species_patterns <- function(species) {
  species <- match.arg(
    species %||% "auto",
    c(
      "auto", "human", "mouse", "rat", "pig", "cow", "cattle",
      "sheep", "goat", "chicken", "duck", "goose", "mixed"
    )
  )
  switch(species,
    human = c("^MT-", "^MT\\.", "^MT_"),
    mouse = c("^mt-", "^Mt-", "^MT-", "^mt\\.", "^MT\\.", "^mt_", "^MT_"),
    rat = c("^mt-", "^Mt-", "^MT-", "^mt\\.", "^MT\\.", "^mt_", "^MT_"),
    pig = common_mito_patterns(),
    cow = common_mito_patterns(),
    cattle = common_mito_patterns(),
    sheep = common_mito_patterns(),
    goat = common_mito_patterns(),
    chicken = common_mito_patterns(),
    duck = common_mito_patterns(),
    goose = common_mito_patterns(),
    mixed = common_mito_patterns(),
    auto = common_mito_patterns()
  )
}

match_exact_features <- function(feature_names, exact_features, ignore_case = FALSE) {
  if (!length(exact_features)) {
    return(character(0))
  }
  if (ignore_case) {
    feature_lookup <- toupper(feature_names)
    exact_lookup <- toupper(exact_features)
    return(unique(feature_names[feature_lookup %in% exact_lookup]))
  }
  unique(intersect(exact_features, feature_names))
}

.validate_existing_mito_ratio <- function(object,
                                          mito_col,
                                          convert_percent = TRUE) {
  values <- object@meta.data[[mito_col]]
  if (!is.numeric(values)) {
    stop("Existing mitochondrial column '", mito_col, "' is not numeric.", call. = FALSE)
  }
  if (any(!is.finite(values))) {
    stop(
      "Existing mitochondrial column '",
      mito_col,
      "' contains NA, NaN, or infinite values.",
      call. = FALSE
    )
  }
  if (length(values) && min(values, na.rm = TRUE) < 0) {
    stop("Existing mitochondrial column '", mito_col, "' contains negative values.", call. = FALSE)
  }

  converted <- FALSE
  max_value <- if (length(values)) max(values, na.rm = TRUE) else 0
  if (is.finite(max_value) && max_value > 1) {
    if (max_value <= 100 && isTRUE(convert_percent)) {
      object@meta.data[[mito_col]] <- values / 100
      converted <- TRUE
    } else if (max_value <= 100) {
      stop(
        "Existing mitochondrial column '",
        mito_col,
        "' appears to be percent-scale. Set convert_percent = TRUE to convert it to a fraction.",
        call. = FALSE
      )
    } else {
      stop(
        "Existing mitochondrial column '",
        mito_col,
        "' is not on a valid fraction or percentage scale.",
        call. = FALSE
      )
    }
  }

  list(object = object, converted = converted)
}

.get_feature_metadata <- function(object, assay) {
  assay_obj <- object[[assay]]
  meta_features <- tryCatch(
    as.data.frame(assay_obj[[]]),
    error = function(e) NULL
  )
  if (!is.data.frame(meta_features) || !nrow(meta_features)) {
    meta_features <- tryCatch(
      as.data.frame(assay_obj@meta.features),
      error = function(e) NULL
    )
  }
  if (!is.data.frame(meta_features) || !nrow(meta_features)) {
    return(NULL)
  }
  meta_features
}

.match_mito_patterns <- function(values,
                                 patterns,
                                 ignore_case = FALSE) {
  unique(unlist(lapply(patterns, function(current_pattern) {
    grep(
      pattern = current_pattern,
      x = values,
      value = TRUE,
      ignore.case = ignore_case
    )
  }), use.names = FALSE))
}

.detect_mito_features <- function(object,
                                  assay,
                                  species = NULL,
                                  mito_features = NULL,
                                  mito_pattern = NULL) {
  feature_names <- rownames(object[[assay]])
  inspected_metadata <- FALSE
  if (!length(feature_names)) {
    return(list(features = character(0), source = "none", inspected_metadata = inspected_metadata))
  }

  if (!is.null(mito_features)) {
    matched <- unique(intersect(as.character(mito_features), feature_names))
    if (!length(matched)) {
      stop(
        "None of the supplied mitochondrial features were found in the Seurat object.",
        call. = FALSE
      )
    }
    return(list(features = matched, source = "user_supplied_features", inspected_metadata = inspected_metadata))
  }

  patterns <- if (!is.null(mito_pattern)) {
    as.character(mito_pattern)
  } else {
    mito_species_patterns(species %||% "auto")
  }
  exact_features <- mito_species_exact_features(species %||% "auto")

  matched <- unique(c(
    .match_mito_patterns(feature_names, patterns, ignore_case = FALSE),
    match_exact_features(feature_names, exact_features, ignore_case = FALSE)
  ))
  if (length(matched)) {
    return(list(features = matched, source = "feature_names", inspected_metadata = inspected_metadata))
  }

  matched <- unique(c(
    .match_mito_patterns(feature_names, patterns, ignore_case = TRUE),
    match_exact_features(feature_names, exact_features, ignore_case = TRUE)
  ))
  if (length(matched)) {
    return(list(features = matched, source = "feature_names_case_insensitive", inspected_metadata = inspected_metadata))
  }

  feature_metadata <- .get_feature_metadata(object, assay)
  inspected_metadata <- TRUE
  if (is.data.frame(feature_metadata) && nrow(feature_metadata)) {
    candidate_columns <- intersect(
      c("gene_name", "gene_short_name", "symbol", "feature_name", "name"),
      colnames(feature_metadata)
    )
    for (column in candidate_columns) {
      values <- as.character(feature_metadata[[column]])
      matched_values <- unique(c(
        .match_mito_patterns(values, patterns, ignore_case = FALSE),
        match_exact_features(values, exact_features, ignore_case = FALSE)
      ))
      if (!length(matched_values)) {
        matched_values <- unique(c(
          .match_mito_patterns(values, patterns, ignore_case = TRUE),
          match_exact_features(values, exact_features, ignore_case = TRUE)
        ))
      }
      if (length(matched_values)) {
        matched_index <- values %in% matched_values
        matched_features <- unique(feature_names[matched_index])
        matched_features <- matched_features[!is.na(matched_features)]
        if (length(matched_features)) {
          return(list(
            features = matched_features,
            source = paste0("feature_metadata:", column),
            inspected_metadata = inspected_metadata
          ))
        }
      }
    }
  }

  list(features = character(0), source = "none", inspected_metadata = inspected_metadata)
}

.store_mito_ratio_info <- function(object,
                                   mito_col,
                                   assay,
                                   species,
                                   n_mito_features,
                                   mito_feature_source,
                                   converted_from_percent,
                                   mito_features_used = character(0)) {
  object@misc$SCdetMito_mito_ratio_info <- list(
    mito_col = mito_col,
    assay = assay,
    species = species %||% NA_character_,
    n_mito_features = n_mito_features,
    mito_feature_source = mito_feature_source,
    scale = "fraction",
    converted_from_percent = isTRUE(converted_from_percent)
  )
  object@misc$SCdetMito <- object@misc$SCdetMito %||% list()
  object@misc$SCdetMito$mito_ratio <- list(
    mito_features_used = mito_features_used,
    n_mito_features = n_mito_features,
    species = species %||% NA_character_,
    pattern = mito_feature_source,
    mito_column = mito_col,
    assay = assay,
    scale = "fraction",
    converted_from_percent = isTRUE(converted_from_percent)
  )
  object
}

#' Ensure mitochondrial ratio metadata are available
#'
#' @description
#' `ensure_mito_ratio()` is the unified mitochondrial fraction preparation
#' function used by SCdetMito workflows. It validates an existing mitochondrial
#' metadata column or calculates one from mitochondrial features when possible.
#'
#' @details
#' Values are stored as fractions on a 0--1 scale. If an existing column appears
#' to be percent-scale and `convert_percent = TRUE`, it is converted to a
#' fraction. When the column is absent, mitochondrial features are detected from
#' user-supplied `mito_features`, a user-supplied `mito_pattern`, species-aware
#' gene-symbol patterns, or feature metadata columns such as `gene_name`,
#' `gene_short_name`, `symbol`, `feature_name`, and `name`.
#'
#' If mitochondrial features were removed during preprocessing, the ratio cannot
#' be reconstructed from the object alone. In that case, provide a valid
#' `mito_features` vector or a precomputed `mito_col` from raw counts that
#' retained mitochondrial features.
#'
#' @param object A Seurat object.
#' @param mito_col Metadata column used for mitochondrial fractions. Defaults to
#'   `"mitoRatio"`.
#' @param species Optional species preset used for mitochondrial feature
#'   detection.
#' @param assay Assay used to identify features and calculate percentages.
#'   Defaults to the active assay.
#' @param mito_features Optional complete character vector of mitochondrial
#'   features retained in the object.
#' @param mito_pattern Optional regular expression used to identify
#'   mitochondrial features.
#' @param overwrite Whether to recalculate `mito_col` when it already exists.
#'   Defaults to `FALSE`.
#' @param convert_percent Whether percent-scale existing values should be
#'   converted to fractions. Defaults to `TRUE`.
#' @param verbose Whether to emit status messages. Defaults to `TRUE`.
#'
#' @return A Seurat object with `object@meta.data[[mito_col]]` stored as a
#'   numeric mitochondrial fraction and summary information in
#'   `object@misc$SCdetMito_mito_ratio_info`.
#' @export
#'
#' @examples
#' seu <- load_demo_pbmc()
#' seu <- ensure_mito_ratio(seu, mito_col = "mitoRatio")
#'
#' mito_features <- grep("^MT-", rownames(seu), value = TRUE)
#' if (length(mito_features)) {
#'   seu <- ensure_mito_ratio(seu, mito_features = mito_features)
#' }
ensure_mito_ratio <- function(object,
                              mito_col = "mitoRatio",
                              species = NULL,
                              assay = NULL,
                              mito_features = NULL,
                              mito_pattern = NULL,
                              overwrite = FALSE,
                              convert_percent = TRUE,
                              verbose = TRUE) {
  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object.", call. = FALSE)
  }
  if (!is.character(mito_col) || length(mito_col) != 1 || !nzchar(mito_col)) {
    stop("'mito_col' must be a single non-empty metadata column name.", call. = FALSE)
  }
  assay <- assay %||% Seurat::DefaultAssay(object)
  if (!assay %in% names(object@assays)) {
    stop("Assay '", assay, "' is not present in the Seurat object.", call. = FALSE)
  }

  if (mito_col %in% colnames(object@meta.data) && !isTRUE(overwrite)) {
    validated <- .validate_existing_mito_ratio(
      object = object,
      mito_col = mito_col,
      convert_percent = convert_percent
    )
    object <- validated$object
    object <- .store_mito_ratio_info(
      object = object,
      mito_col = mito_col,
      assay = assay,
      species = species,
      n_mito_features = NA_integer_,
      mito_feature_source = "existing_metadata_column",
      converted_from_percent = validated$converted,
      mito_features_used = character(0)
    )
    if (isTRUE(verbose)) {
      message("Using existing mitochondrial fraction column: ", mito_col)
    }
    return(object)
  }

  detection <- .detect_mito_features(
    object = object,
    assay = assay,
    species = species,
    mito_features = mito_features,
    mito_pattern = mito_pattern
  )
  if (!length(detection$features)) {
    feature_names <- rownames(object[[assay]])
    stop(
      "No mitochondrial features were detected. Assay: ",
      assay,
      ". Number of features: ",
      length(feature_names),
      ". First 20 feature names: ",
      paste(utils::head(feature_names, 20), collapse = ", "),
      ". Species: ",
      species %||% "not supplied",
      ". Feature metadata inspected: ",
      if (isTRUE(detection$inspected_metadata)) "yes" else "no",
      ". Please provide mitochondrial features manually with 'mito_features' or use raw counts retaining mitochondrial genes.",
      call. = FALSE
    )
  }

  percent_values <- Seurat::PercentageFeatureSet(
    object = object,
    features = detection$features,
    assay = assay
  )
  object@meta.data[[mito_col]] <- percent_values / 100
  object <- .store_mito_ratio_info(
    object = object,
    mito_col = mito_col,
    assay = assay,
    species = species,
    n_mito_features = length(detection$features),
    mito_feature_source = detection$source,
    converted_from_percent = FALSE,
    mito_features_used = detection$features
  )
  if (isTRUE(verbose)) {
    message(
      "Calculated mitochondrial fraction column '",
      mito_col,
      "' using ",
      length(detection$features),
      " mitochondrial features."
    )
  }
  object
}

#' add_mitoRatio: backward-compatible mitochondrial ratio wrapper
#'
#' @description
#' `add_mitoRatio()` is retained for backward compatibility. New workflows can
#' use [ensure_mito_ratio()] directly.
#'
#' @details
#' The function delegates mitochondrial feature detection and fraction
#' calculation to [ensure_mito_ratio()]. Manual mitochondrial feature lists
#' should include all mitochondrial features retained in the object. Short
#' illustrative subsets should not be used for real QC.
#'
#' @param seurat_obj A Seurat object.
#' @param pattern Regex pattern used to match mitochondrial features. Use
#'   `"auto"` to scan common mitochondrial prefixes. Defaults to `"auto"`.
#' @param features Optional complete character vector of mitochondrial feature
#'   names. When supplied, `pattern` and `species` are ignored.
#' @param species Optional species preset used for mitochondrial feature
#'   detection.
#' @param assay Assay used to calculate feature percentages. Defaults to the
#'   active assay.
#' @param column Metadata column name used to store the result. Defaults to
#'   `"mitoRatio"`.
#' @param scale Output scale for the stored metric. `"fraction"` stores values
#'   between 0 and 1 and `"percent"` stores values between 0 and 100.
#' @param ignore_case Retained for backward compatibility. Matching now tries
#'   case-sensitive and case-insensitive detection internally.
#'
#' @return The input Seurat object with an added mitochondrial ratio column in
#'   `meta.data`.
#' @export
add_mitoRatio <- function(seurat_obj,
                          pattern = "auto",
                          features = NULL,
                          species = NULL,
                          assay = NULL,
                          column = "mitoRatio",
                          scale = c("fraction", "percent"),
                          ignore_case = FALSE) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.", call. = FALSE)
  }
  scale <- match.arg(scale)
  mito_pattern <- if (identical(pattern, "auto")) NULL else pattern
  updated <- ensure_mito_ratio(
    object = seurat_obj,
    mito_col = column,
    species = species,
    assay = assay,
    mito_features = features,
    mito_pattern = mito_pattern,
    overwrite = TRUE,
    convert_percent = TRUE,
    verbose = FALSE
  )
  if (identical(scale, "percent")) {
    updated@meta.data[[column]] <- updated@meta.data[[column]] * 100
    updated@misc$SCdetMito_mito_ratio_info$scale <- "percent"
    updated@misc$SCdetMito$mito_ratio$scale <- "percent"
  }
  updated
}
