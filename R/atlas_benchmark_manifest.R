#' atlas_benchmark_manifest: curated benchmark datasets for broader validation
#'
#' @description
#' `atlas_benchmark_manifest()` returns a curated table of public scRNA-seq or
#' snRNA-seq datasets that can be considered for validating `SCdetMito` across
#' authoritative lung cohorts, multi-organ atlases, and cross-species studies.
#'
#' @details
#' The manifest is intended to support manuscript planning and reproducible
#' validation design. It distinguishes between:
#'
#' - `main_text`: strong candidates for the primary paper benchmark
#' - `independent_validation`: additional cohorts for replication
#' - `supplementary_extension`: broader atlas studies for cross-species or
#'   multi-tissue generalization
#'
#' This is a reference manifest only. Returning a row does not validate current
#' dataset availability, download completeness, annotation quality, or local
#' reusability. Users must independently verify accessions, license terms,
#' data integrity, and analysis-ready input files before using any dataset for
#' manuscript claims.
#'
#' `accession_or_portal` may be a GEO/BioProject accession, a portal name, or
#' `NA` when the source is known but the exact packaged input path is left to a
#' local curation step.
#'
#' @param role Optional role filter. Defaults to `NULL`, which returns all rows.
#' @param species Optional species filter. Defaults to `NULL`.
#' @param tissue Optional tissue filter. Defaults to `NULL`.
#'
#' @return A data frame describing recommended benchmark candidates.
#' @export
atlas_benchmark_manifest <- function(role = NULL, species = NULL, tissue = NULL) {
  manifest <- data.frame(
    dataset_id = c(
      "GSE135893",
      "GSE136831",
      "GSE157996",
      "GSE159354",
      "GSE159929",
      "GSE200090",
      "GSE307046",
      "GSE233285",
      "GSE176512",
      "GSE158117",
      "GSE184343",
      "GSE224329",
      "GSE183300"
    ),
    species = c(
      "human", "human", "human", "human", "human",
      "mouse", "rat", "pig", "cattle", "cattle", "sheep", "chicken", "duck"
    ),
    tissue = c(
      "lung",
      "lung",
      "lung",
      "lung",
      "multi-organ",
      "multi-organ",
      "kidney",
      "multi-organ",
      "liver",
      "mammary",
      "testis",
      "PBMC",
      "lung"
    ),
    disease_context = c(
      "IPF vs control",
      "IPF/COPD/control",
      "IPF vs healthy",
      "ILD/IPF/control",
      "healthy reference",
      "healthy reference",
      "CIH/OSA vs normoxia",
      "healthy reference",
      "healthy liver reference",
      "development/lactation reference",
      "development reference",
      "healthy immune reference",
      "cross-species lung reference"
    ),
    journal = c(
      "Science Advances",
      "Science Advances",
      "Journal of Clinical Investigation",
      "mixed downstream publications",
      "Genome Biology",
      "Communications Biology",
      "Scientific Reports",
      "bioRxiv / GEO atlas",
      "GEO atlas",
      "citation missing / GEO atlas",
      "Current Issues in Molecular Biology",
      "BMC Genomics",
      "Nature Communications"
    ),
    year = c(2020, 2020, 2022, 2020, 2020, 2022, 2025, 2025, 2022, 2023, 2022, 2024, 2021),
    accession_or_portal = c(
      "GSE135893",
      "GSE136831",
      "GSE157996",
      "GSE159354",
      "GSE159929",
      "GSE200090",
      "GSE307046",
      "GSE233285",
      "GSE176512",
      "GSE158117",
      "GSE184343",
      "GSE224329",
      "GSE183300"
    ),
    role = rep(c("main_text", rep("independent_validation", 3), rep("supplementary_extension", 9)), length.out = 13),
    has_raw_counts = rep(TRUE, 13),
    has_author_annotations = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    availability_validated = rep(FALSE, 13),
    manuscript_priority = c(
      "highest",
      "high",
      "medium",
      "medium",
      "high",
      "medium",
      "exploratory",
      "medium",
      "exploratory",
      "exploratory",
      "exploratory",
      "exploratory",
      "medium"
    ),
    status = c(
      "in_progress",
      "planned",
      "executed",
      "planned",
      "executed",
      "executed",
      "executed",
      "executed",
      "executed",
      "executed",
      "executed",
      "in_progress",
      "in_progress"
    ),
    notes = c(
      "Strong candidate for a direct author-annotation versus SCdetMito QC benchmark in human IPF lung.",
      "Large authoritative lung atlas suitable as a second major replication cohort.",
      "Already operational in the repository; useful as a lightweight real-data application set.",
      "Useful epithelial validation cohort, but less persuasive as the single flagship benchmark.",
      "Adult human cell atlas spanning 15 organs with public count matrices for multi-tissue transfer validation.",
      "Mouse atlas with lung, kidney, and spleen coverage suitable for rodent transfer checks.",
      "Rat kidney single-nucleus dataset useful for confirming species transfer beyond mouse.",
      "Pig multi-tissue atlas covering adipose, heart, kidney, liver, muscle, and spleen for livestock validation.",
      "Bovine liver route with public MTX triplets; now executed as a strict target-tissue cattle validation.",
      "Bovine mammary scRNA-seq atlas offering a pragmatic large-animal extension with public MTX files.",
      "Sheep testis scRNA-seq atlas with reusable MTX files; now executed as an extra ruminant transfer route in the repository.",
      "Chicken PBMC dataset with public matrix files; useful to test avian mitochondrial feature handling.",
      "Nature Communications cross-species atlas including duck lung, useful for avian lung transfer validation."
    ),
    stringsAsFactors = FALSE
  )

  if (!is.null(role)) {
    manifest <- manifest[manifest$role %in% role, , drop = FALSE]
  }
  if (!is.null(species)) {
    manifest <- manifest[manifest$species %in% species, , drop = FALSE]
  }
  if (!is.null(tissue)) {
    manifest <- manifest[manifest$tissue %in% tissue, , drop = FALSE]
  }

  rownames(manifest) <- NULL
  manifest
}
