common_mito_patterns <- function() {
  c("^MT-", "^mt-", "^Mt-", "^mt:")
}

mito_species_exact_features <- function(species) {
  species <- match.arg(
    species,
    c(
      "auto", "human", "mouse", "rat", "pig", "cow", "cattle",
      "sheep", "goat", "chicken", "duck", "goose", "mixed"
    )
  )

  bare_mt_symbols <- c(
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6",
    "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"
  )

  switch(species,
    human = c(
      "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8",
      "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4",
      "MT-ND5", "MT-ND6", "MT-CYB"
    ),
    mouse = c(
      "mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp8",
      "mt-Atp6", "mt-Co3", "mt-Nd3", "mt-Nd4l", "mt-Nd4",
      "mt-Nd5", "mt-Nd6", "mt-Cytb",
      "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8",
      "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4",
      "MT-ND5", "MT-ND6", "MT-CYB"
    ),
    rat = c(
      "mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp8",
      "mt-Atp6", "mt-Co3", "mt-Nd3", "mt-Nd4l", "mt-Nd4",
      "mt-Nd5", "mt-Nd6", "mt-Cytb",
      "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8",
      "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4",
      "MT-ND5", "MT-ND6", "MT-CYB"
    ),
    pig = bare_mt_symbols,
    cow = bare_mt_symbols,
    cattle = bare_mt_symbols,
    sheep = bare_mt_symbols,
    goat = bare_mt_symbols,
    chicken = bare_mt_symbols,
    duck = bare_mt_symbols,
    goose = bare_mt_symbols,
    mixed = unique(c(
      bare_mt_symbols,
      "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8",
      "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4",
      "MT-ND5", "MT-ND6", "MT-CYB",
      "mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp8",
      "mt-Atp6", "mt-Co3", "mt-Nd3", "mt-Nd4l", "mt-Nd4",
      "mt-Nd5", "mt-Nd6", "mt-Cytb"
    )),
    auto = character(0)
  )
}

mito_species_patterns <- function(species) {
  species <- match.arg(
    species,
    c(
      "auto", "human", "mouse", "rat", "pig", "cow", "cattle",
      "sheep", "goat", "chicken", "duck", "goose", "mixed"
    )
  )
  switch(species,
    human = "^MT-",
    mouse = c("^mt-", "^Mt-", "^MT-"),
    rat = c("^mt-", "^Mt-", "^MT-"),
    pig = c("^MT-", "^mt-", "^Mt-"),
    cow = c("^MT-", "^mt-", "^Mt-"),
    cattle = c("^MT-", "^mt-", "^Mt-"),
    sheep = c("^MT-", "^mt-", "^Mt-"),
    goat = c("^MT-", "^mt-", "^Mt-"),
    chicken = c("^MT-", "^mt-", "^Mt-"),
    duck = c("^MT-", "^mt-", "^Mt-"),
    goose = c("^MT-", "^mt-", "^Mt-"),
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

  intersect(exact_features, feature_names)
}

resolve_mito_features <- function(feature_names,
                                  pattern = "auto",
                                  features = NULL,
                                  species = NULL,
                                  ignore_case = FALSE) {
  if (!is.null(features)) {
    matched_features <- intersect(as.character(features), feature_names)
    if (!length(matched_features)) {
      stop(
        "None of the requested mitochondrial features were found in the Seurat object.",
        call. = FALSE
      )
    }
    return(unique(matched_features))
  }

  candidate_patterns <- if (!is.null(species)) {
    mito_species_patterns(species)
  } else if (identical(pattern, "auto")) {
    common_mito_patterns()
  } else {
    as.character(pattern)
  }
  exact_features <- if (!is.null(species)) {
    mito_species_exact_features(species)
  } else {
    character(0)
  }

  matched_by_pattern <- unique(unlist(lapply(candidate_patterns, function(current_pattern) {
    grep(
      pattern = current_pattern,
      x = feature_names,
      value = TRUE,
      ignore.case = ignore_case
    )
  }), use.names = FALSE))
  matched_by_exact <- match_exact_features(
    feature_names = feature_names,
    exact_features = exact_features,
    ignore_case = ignore_case
  )
  matched_features <- unique(c(matched_by_pattern, matched_by_exact))

  if (!length(matched_features)) {
    stop(
      "No mitochondrial features were detected. Please check gene naming style, species, or provide mitochondrial features manually.",
      call. = FALSE
    )
  }

  matched_features
}

#' add_mitoRatio: calculate and store mitochondrial ratios in a Seurat object
#'
#' @description
#' `add_mitoRatio()` calculates mitochondrial abundance for each cell and stores
#' the result in `meta.data`.
#'
#' @details
#' By default the function uses `pattern = "auto"` and searches for common
#' mitochondrial prefixes such as `MT-`, `mt-`, and `Mt-`. This makes the
#' function safer for public data reuse and cross-species workflows where gene
#' naming conventions vary. When `species` is supplied, `add_mitoRatio()`
#' additionally tries curated exact mitochondrial gene symbols for that species.
#' This is important for livestock datasets where mitochondrial genes are often
#' annotated as symbols such as `ND1`, `COX1`, `ATP6`, and `CYTB` rather than
#' `MT-`-prefixed features.
#'
#' @param seurat_obj A Seurat object.
#' @param pattern Regex pattern used to match mitochondrial features. Use
#'   `"auto"` to scan common mitochondrial prefixes. Defaults to `"auto"`.
#' @param features Optional character vector of mitochondrial feature names. When
#'   supplied, `pattern` and `species` are ignored.
#' @param species Optional species preset used for mitochondrial feature
#'   detection. Supported values are `"auto"`, `"human"`, `"mouse"`, `"rat"`,
#'   `"pig"`, `"cow"`, `"cattle"`, `"sheep"`, `"goat"`, `"chicken"`,
#'   `"duck"`, `"goose"`, and `"mixed"`. Defaults to `NULL`.
#' @param assay Assay used to calculate feature percentages. Defaults to the
#'   active assay.
#' @param column Metadata column name used to store the result. Defaults to
#'   `"mitoRatio"`.
#' @param scale Output scale for the stored metric. `"fraction"` stores values
#'   between 0 and 1 and `"percent"` stores values between 0 and 100.
#' @param ignore_case Whether feature matching should ignore case when `pattern`
#'   or `species`-based detection is used. Defaults to `FALSE`.
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
  assay <- if (is.null(assay)) Seurat::DefaultAssay(seurat_obj) else assay
  feature_names <- rownames(seurat_obj[[assay]])
  mito_features <- resolve_mito_features(
    feature_names = feature_names,
    pattern = pattern,
    features = features,
    species = species,
    ignore_case = ignore_case
  )

  percent_values <- Seurat::PercentageFeatureSet(
    object = seurat_obj,
    features = mito_features,
    assay = assay
  )

  seurat_obj@meta.data[[column]] <- if (identical(scale, "fraction")) {
    percent_values / 100
  } else {
    percent_values
  }

  seurat_obj@misc$SCdetMito <- seurat_obj@misc$SCdetMito %||% list()
  seurat_obj@misc$SCdetMito$mito_ratio <- list(
    mito_features_used = mito_features,
    n_mito_features = length(mito_features),
    species = species %||% NA_character_,
    pattern = if (length(pattern)) paste(pattern, collapse = ";") else NA_character_,
    mito_column = column,
    assay = assay,
    scale = scale
  )

  seurat_obj
}
