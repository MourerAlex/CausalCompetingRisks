#' Compute Causal Contrasts
#'
#' Computes risk differences and risk ratios for total, separable direct,
#' and separable indirect effects from cumulative incidence estimates.
#'
#' @param cumulative_incidence Data frame with columns `k`, `arm_11`,
#'   `arm_00`, `arm_10`.
#' @param eval_times Numeric vector or NULL. Time points at which to evaluate
#'   contrasts. NULL = all time points.
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{k}{Time point.}
#'     \item{total_rd}{Total effect risk difference: arm(1,1) - arm(0,0).}
#'     \item{total_rr}{Total effect risk ratio: arm(1,1) / arm(0,0).}
#'     \item{sep_direct_rd}{Separable direct (A_Y) RD: arm(1,0) - arm(0,0).}
#'     \item{sep_direct_rr}{Separable direct (A_Y) RR: arm(1,0) / arm(0,0).}
#'     \item{sep_indirect_rd}{Separable indirect (A_D) RD: arm(1,1) - arm(1,0).}
#'     \item{sep_indirect_rr}{Separable indirect (A_D) RR: arm(1,1) / arm(1,0).}
#'   }
#'
#' @family internal
#' @keywords internal
compute_contrasts <- function(cumulative_incidence, eval_times = NULL) {
  ci <- cumulative_incidence

  if (!is.null(eval_times)) {
    ci <- ci[ci$k %in% eval_times, ]
  }

  data.frame(
    k = ci$k,
    # Total effect: (1,1) vs (0,0)
    total_rd = ci$arm_11 - ci$arm_00,
    total_rr = ci$arm_11 / ci$arm_00,
    # Separable direct (A_Y): (1,0) vs (0,0)
    sep_direct_rd = ci$arm_10 - ci$arm_00,
    sep_direct_rr = ci$arm_10 / ci$arm_00,
    # Separable indirect (A_D): (1,1) vs (1,0)
    sep_indirect_rd = ci$arm_11 - ci$arm_10,
    sep_indirect_rr = ci$arm_11 / ci$arm_10,
    row.names = NULL
  )
}


#' Build Cumulative Incidence from Results
#'
#' Combines g-formula and/or IPW cumulative incidence results into a single
#' data frame.
#'
#' @param results List with `$gformula` and/or `$ipw` elements.
#' @param method Character. Estimation method used.
#'
#' @return Data frame with cumulative incidence.
#'
#' @family internal
#' @keywords internal
build_cumulative_incidence <- function(results, method) {
  if (method == "gformula") {
    results$gformula
  } else if (method %in% c("ipw1", "ipw2")) {
    results$ipw$cumulative_incidence
  } else {
    # method = "both": return g-formula as primary, IPW stored in diagnostics
    results$gformula
  }
}


#' Build Diagnostics from Results
#'
#' Extracts diagnostic information (weight summaries, truncated IDs) from
#' IPW results.
#'
#' @param results List with estimation results.
#' @param method Character. Estimation method used.
#'
#' @return List with diagnostic information, or NULL.
#'
#' @family internal
#' @keywords internal
build_diagnostics <- function(results, method) {
  if (method %in% c("both", "ipw1", "ipw2") && !is.null(results$ipw)) {
    list(
      weight_summary = results$ipw$weight_summary,
      truncated_ids = results$ipw$truncated_ids,
      ipw_cumulative_incidence = results$ipw$cumulative_incidence
    )
  } else {
    NULL
  }
}
