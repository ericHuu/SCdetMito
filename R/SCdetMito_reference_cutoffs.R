# SCdetMito
# Author: Silu Hu
# Contact: husilu0902@gmail.com
# Version: 1.4.1
# Last updated: 2026-05-21

#' Literature-informed mitochondrial cutoff references
#'
#' @description
#' `SCdetMito_reference_cutoffs()` returns a lightweight table of
#' literature-informed mitochondrial fraction cutoffs used as interpretive
#' references by SCdetMito. These values are QC decision priors, not hard rules
#' or universal biological thresholds.
#'
#' @details
#' The table is intentionally conservative. Species-level global references are
#' provided where tissue-specific numeric values should not be over-specified.
#' Tissue rows add context notes for interpretation, especially when tissue
#' handling, dissociation stress, post-mortem interval, disease state, or cell
#' composition can shift mitochondrial distributions.
#'
#' @param species Optional species filter. Matching is case-insensitive.
#' @param tissue Optional tissue filter. Matching is case-insensitive. If a
#'   species is supplied and no exact tissue match exists, the species-level
#'   global row is returned with an informational message.
#'
#' @return A data.frame with columns `species`, `tissue`, `cell_context`,
#'   `reference_cutoff`, `unit`, `evidence_level`, `reference_scope`,
#'   `is_tissue_specific`, `source_short`, `source_detail`, `doi_or_url`, and
#'   `note`.
#' @export
#'
#' @examples
#' SCdetMito_reference_cutoffs()
#' SCdetMito_reference_cutoffs(species = "human")
#' SCdetMito_reference_cutoffs(species = "mouse", tissue = "heart")
SCdetMito_reference_cutoffs <- function(species = NULL, tissue = NULL) {
  ref <- data.frame(
    species = c("human", "mouse", "human", "mouse", "mouse", "human"),
    tissue = c("global", "global", "adipose", "heart", "kidney", "PBMC"),
    cell_context = c(
      "scRNA-seq",
      "scRNA-seq",
      "adipose tissue scRNA-seq",
      "heart tissue scRNA-seq",
      "kidney tissue scRNA-seq",
      "peripheral blood mononuclear cells"
    ),
    reference_cutoff = c(0.10, 0.05, 0.10, 0.05, 0.05, 0.10),
    unit = rep("fraction", 6),
    evidence_level = c(
      "species-level literature prior",
      "species-level literature prior",
      "species-level prior with tissue caution",
      "species-level prior with tissue caution",
      "species-level prior with tissue caution",
      "species-level prior with PBMC caution"
    ),
    reference_scope = c(
      "species_global_prior",
      "species_global_prior",
      "species_global_prior_with_tissue_note",
      "species_global_prior_with_tissue_note",
      "species_global_prior_with_tissue_note",
      "species_global_prior_with_tissue_note"
    ),
    is_tissue_specific = rep(FALSE, 6),
    source_short = c(
      "Osorio and Cai 2021",
      "Osorio and Cai 2021",
      "Osorio and Cai 2021",
      "Osorio and Cai 2021",
      "Osorio and Cai 2021",
      "Osorio and Cai 2021"
    ),
    source_detail = c(
      "Osorio and Cai, Systematic determination of the mitochondrial proportion in human and mice tissues for single-cell RNA-sequencing data quality control, Bioinformatics, 2021.",
      "Osorio and Cai, Systematic determination of the mitochondrial proportion in human and mice tissues for single-cell RNA-sequencing data quality control, Bioinformatics, 2021.",
      "Osorio and Cai, Bioinformatics, 2021; adipose and other metabolically active tissues may exceed simple global assumptions and require inspection.",
      "Osorio and Cai, Bioinformatics, 2021; heart can show high mitochondrial background and data-driven boundaries may deviate from global references.",
      "Osorio and Cai, Bioinformatics, 2021; kidney mitochondrial distributions can be tissue- and protocol-dependent.",
      "Osorio and Cai, Bioinformatics, 2021; PBMCs often have relatively low mitochondrial fractions, but inspection remains required."
    ),
    doi_or_url = rep("https://doi.org/10.1093/bioinformatics/btaa751", 6),
    note = c(
      "Human global reference used as a decision prior; not a hard filtering rule.",
      "Mouse global reference used as a decision prior; not a hard filtering rule.",
      "This row applies a species-level prior with tissue-specific interpretation notes; it is not an exact tissue-specific cutoff. Inspect adipose tissue distributions before filtering.",
      "This row applies a species-level prior with tissue-specific interpretation notes; it is not an exact tissue-specific cutoff. High-mitochondrial burden can reflect tissue biology, handling, or cell stress.",
      "This row applies a species-level prior with tissue-specific interpretation notes; it is not an exact tissue-specific cutoff. Inspect retained-cell and interval-loss profiles.",
      "This row applies a species-level prior with tissue-specific interpretation notes; it is not an exact tissue-specific cutoff. PBMC references are dataset-specific; 5% may be stringent in some workflows and 10% is retained here as a cautious upper prior."
    ),
    stringsAsFactors = FALSE
  )

  out <- ref
  if (!is.null(species)) {
    species_key <- tolower(as.character(species[[1]]))
    out <- out[tolower(out$species) == species_key, , drop = FALSE]
  }
  if (!is.null(tissue)) {
    tissue_key <- tolower(as.character(tissue[[1]]))
    exact <- out[tolower(out$tissue) == tissue_key, , drop = FALSE]
    if (nrow(exact)) {
      out <- exact
    } else if (!is.null(species)) {
      global <- out[tolower(out$tissue) == "global", , drop = FALSE]
      if (nrow(global)) {
        message("No exact tissue-specific reference found; returning species-level global reference.")
        out <- global
      } else {
        out <- exact
      }
    } else {
      out <- exact
    }
  }

  rownames(out) <- NULL
  out
}

.get_reference_cutoff <- function(species = NULL,
                                  tissue = NULL,
                                  default_reference = NULL,
                                  reference_table = NULL,
                                  verbose = TRUE) {
  if (!is.null(default_reference)) {
    reference_value <- normalize_mito_cutoff_value(default_reference, "reference_cutoff")
    return(list(
      reference_cutoff = reference_value,
      reference_species = if (is.null(species)) NA_character_ else as.character(species[[1]]),
      reference_tissue = if (is.null(tissue)) NA_character_ else as.character(tissue[[1]]),
      reference_source = "user_provided",
      reference_note = "User-provided literature/reference cutoff.",
      reference_found = TRUE
    ))
  }

  ref <- if (is.null(reference_table)) {
    SCdetMito_reference_cutoffs()
  } else {
    reference_table
  }
  required <- c("species", "tissue", "reference_cutoff")
  if (!all(required %in% colnames(ref))) {
    stop(
      "'reference_table' must contain columns: ",
      paste(required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (is.null(species) && is.null(tissue)) {
    return(list(
      reference_cutoff = NA_real_,
      reference_species = NA_character_,
      reference_tissue = NA_character_,
      reference_source = NA_character_,
      reference_note = "No species or tissue reference requested.",
      reference_found = FALSE
    ))
  }

  candidates <- ref
  if (!is.null(species)) {
    species_key <- tolower(as.character(species[[1]]))
    candidates <- candidates[tolower(candidates$species) == species_key, , drop = FALSE]
  }
  exact <- candidates
  if (!is.null(tissue)) {
    tissue_key <- tolower(as.character(tissue[[1]]))
    exact <- candidates[tolower(candidates$tissue) == tissue_key, , drop = FALSE]
  }

  selected <- exact
  if (!nrow(selected) && !is.null(species)) {
    selected <- candidates[tolower(candidates$tissue) == "global", , drop = FALSE]
    if (nrow(selected) && isTRUE(verbose) && !is.null(tissue)) {
      message("No exact tissue-specific reference found; returning species-level global reference.")
    }
  }

  if (!nrow(selected)) {
    return(list(
      reference_cutoff = NA_real_,
      reference_species = if (is.null(species)) NA_character_ else as.character(species[[1]]),
      reference_tissue = if (is.null(tissue)) NA_character_ else as.character(tissue[[1]]),
      reference_source = NA_character_,
      reference_note = "No matching reference cutoff found.",
      reference_found = FALSE
    ))
  }

  selected <- selected[1, , drop = FALSE]
  list(
    reference_cutoff = normalize_mito_cutoff_value(selected$reference_cutoff[[1]], "reference_cutoff"),
    reference_species = selected$species[[1]],
    reference_tissue = selected$tissue[[1]],
    reference_source = if ("source_short" %in% colnames(selected)) selected$source_short[[1]] else NA_character_,
    reference_note = if ("note" %in% colnames(selected)) selected$note[[1]] else NA_character_,
    reference_found = TRUE
  )
}
