#' SCdetMito_methods: list available cutoff-detection options
#'
#' @description
#' `SCdetMito_methods()` returns the currently supported statistical tests,
#' multiple-testing correction methods, and sample-level cutoff selection rules
#' used by `SCdetMito()`.
#'
#' @details
#' The default loss test is `mad_zscore`, a robust upper-tail detector based on
#' the median and MAD of the same-sample interval-specific cell losses. It is
#' intended as the primary FDR-compatible mode. `empirical_tail` remains
#' available as a non-parametric, highly discrete sensitivity/exploratory mode;
#' `poisson_tail` and `zscore` are parametric sensitivity checks, and
#' `threshold_only` bypasses statistical enrichment for explicit rule-based QC.
#'
#' @return A named list of data frames describing selectable methods.
#' @export
SCdetMito_methods <- function() {
  list(
    loss_tests = data.frame(
      method = c("mad_zscore", "empirical_tail", "poisson_tail", "zscore", "threshold_only"),
      default = c(TRUE, FALSE, FALSE, FALSE, FALSE),
      description = c(
        "Default robust normal upper-tail test using median and median absolute deviation of same-sample interval-specific losses.",
        "Non-parametric rank-based upper-tail probability from observed background interval losses; highly discrete under dense scanning.",
        "Poisson upper-tail test against the mean background interval loss.",
        "Normal upper-tail test using background-loss mean and standard deviation.",
        "No enrichment model; treats any qualified loss as significant and is mainly for sensitivity checks."
      ),
      stringsAsFactors = FALSE
    ),
    p_adjust_methods = data.frame(
      method = stats::p.adjust.methods,
      default = stats::p.adjust.methods == "BH",
      stringsAsFactors = FALSE
    ),
    sample_cutoff_methods = data.frame(
      method = c("largest_drop", "max_significant", "median_significant", "min_significant"),
      default = c(TRUE, FALSE, FALSE, FALSE),
      description = c(
        "FDR-significant cutoff with the largest retained-cell loss interval.",
        "Largest FDR-significant cutoff for each sample.",
        "Median FDR-significant cutoff snapped to the cutoff grid.",
        "Smallest FDR-significant cutoff for each sample."
      ),
      stringsAsFactors = FALSE
    ),
    final_cutoff_rule = data.frame(
      rule = "sample_supported_global_cutoff",
      default = TRUE,
      description = "Return the largest global cutoff supported by at least half of detected sample-level cutoffs by default.",
      stringsAsFactors = FALSE
    )
  )
}
