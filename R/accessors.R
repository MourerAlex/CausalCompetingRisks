#' Extract Cumulative Incidence Curves
#'
#' Returns a `"causal_cr_risk"` object containing cumulative incidence
#' estimates for each treatment arm. Has its own `print()` and `plot()`
#' methods.
#'
#' @param fit A `"causal_cr"` object from [causal_cr()].
#'
#' @return An S3 object of class `"causal_cr_risk"` with:
#'   \describe{
#'     \item{cumulative_incidence}{Data frame with `k`, `arm_11`, `arm_00`,
#'       `arm_10`.}
#'     \item{bootstrap}{CI bands if bootstrap was used, else NULL.}
#'   }
#'
#' @seealso [contrast()], [diagnostic()], [plot.causal_cr_risk()]
#' @family accessors
#' @export
risk <- function(fit) {
  stopifnot(inherits(fit, "causal_cr"))
  structure(
    list(
      cumulative_incidence = fit$cumulative_incidence,
      bootstrap = if (!is.null(fit$bootstrap)) fit$bootstrap$ci_curves else NULL
    ),
    class = "causal_cr_risk"
  )
}


#' Extract Causal Contrasts
#'
#' Returns a `"causal_cr_contrast"` object containing risk differences and
#' risk ratios for total, separable direct, and separable indirect effects.
#' Has its own `print()` and `plot()` methods.
#'
#' @param fit A `"causal_cr"` object from [causal_cr()].
#'
#' @return An S3 object of class `"causal_cr_contrast"` with:
#'   \describe{
#'     \item{contrasts}{Data frame with RD and RR at each time point.}
#'     \item{bootstrap}{CI on contrasts if bootstrap was used, else NULL.}
#'   }
#'
#' @seealso [risk()], [diagnostic()], [plot.causal_cr_contrast()]
#' @family accessors
#' @export
contrast <- function(fit) {
  stopifnot(inherits(fit, "causal_cr"))
  structure(
    list(
      contrasts = fit$contrasts,
      bootstrap = if (!is.null(fit$bootstrap)) fit$bootstrap$ci_contrasts else NULL
    ),
    class = "causal_cr_contrast"
  )
}


#' Extract Diagnostics
#'
#' Returns a `"causal_cr_diagnostic"` object containing weight summaries,
#' isolation checks, and flagged IDs. Has its own `print()` and `plot()`
#' methods.
#'
#' @param fit A `"causal_cr"` object from [causal_cr()].
#'
#' @return An S3 object of class `"causal_cr_diagnostic"` with:
#'   \describe{
#'     \item{weight_summary}{Data frame summarizing weight distributions.}
#'     \item{truncated_ids}{IDs with adjusted weights.}
#'     \item{ipw_cumulative_incidence}{IPW results (when method = "both").}
#'   }
#'
#' @seealso [risk()], [contrast()], [plot.causal_cr_diagnostic()]
#' @family accessors
#' @export
diagnostic <- function(fit) {
  stopifnot(inherits(fit, "causal_cr"))
  structure(
    fit$diagnostics %||% list(
      weight_summary = NULL,
      truncated_ids = integer(0),
      ipw_cumulative_incidence = NULL
    ),
    class = "causal_cr_diagnostic"
  )
}
