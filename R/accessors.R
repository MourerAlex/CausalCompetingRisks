#' Extract Cumulative Incidence Curves
#'
#' Returns a `"separable_effects_risk"` object containing cumulative incidence
#' estimates for each treatment arm, across all methods that were run.
#' Optionally pairs with a bootstrap object to surface confidence bands
#' downstream (in `plot.separable_effects_risk`).
#'
#' @param fit A `"separable_effects"` object from [separable_effects()].
#' @param ci Optional. A `"causal_cr_bootstrap"` object from [bootstrap()].
#'   When provided, its `ci_curves` are attached for plotting.
#'
#' @return An S3 object of class `"separable_effects_risk"` with:
#'   \describe{
#'     \item{cumulative_incidence}{Named list per method of cumulative
#'       incidence data.frames (from `fit$cumulative_incidence`).}
#'     \item{ci_curves}{Named list per method of CI-band data.frames from
#'       the bootstrap, or NULL if `ci` was not supplied.}
#'   }
#'
#' @seealso [contrast()], [diagnostic()], [plot.separable_effects_risk()],
#'   [bootstrap()]
#' @family accessors
#' @export
risk <- function(fit, ci = NULL) {
  stopifnot(inherits(fit, "separable_effects"))
  if (!is.null(ci)) {
    stopifnot(inherits(ci, "causal_cr_bootstrap"))
  }
  structure(
    list(
      cumulative_incidence = fit$cumulative_incidence,
      ci_curves            = if (!is.null(ci)) ci$ci_curves else NULL,
      # Full bootstrap replicates + alpha for PROPER per-contrast CIs in
      # plot.separable_effects_risk (per-replicate difference, then quantile).
      replicates           = if (!is.null(ci)) ci$replicates else NULL,
      alpha                = if (!is.null(ci)) ci$alpha      else NULL,
      # References needed for plot's risk_table option (avoids forcing
      # the user to pass `fit` again at plot time). These are R pointers,
      # not copies — no real memory overhead unless the underlying data
      # is modified.
      person_time          = fit$person_time,
      id_col               = fit$id_col,
      treatment_col        = fit$treatment_col,
      times                = fit$times
    ),
    class = "separable_effects_risk"
  )
}


#' Extract Causal Contrasts (Long Format with Confidence Intervals)
#'
#' Returns a `"causal_cr_contrast"` object with a long-format data frame
#' of total, separable direct, and separable indirect effects for a single
#' estimation method. Bootstrap confidence intervals are required — the
#' package does not report contrasts without uncertainty.
#'
#' @param fit A `"separable_effects"` object from [separable_effects()].
#' @param method Character (length 1). Which method's cumulative incidence
#'   to contrast. Must be in `names(fit$cumulative_incidence)`.
#' @param ci A `"causal_cr_bootstrap"` object from [bootstrap()]. Required.
#'
#' @return An S3 object of class `"causal_cr_contrast"` with:
#'   \describe{
#'     \item{contrasts}{Long-format data frame with columns `k`, `contrast`,
#'       `decomp`, `measure`, `estimate`, `lower`, `upper`. 10 rows per
#'       time point (1 total + 2 direct + 2 indirect, RD and RR each).}
#'     \item{method}{The method used.}
#'     \item{alpha}{Significance level from the bootstrap object.}
#'   }
#'
#' @seealso [risk()], [diagnostic()], [plot.causal_cr_contrast()],
#'   [bootstrap()]
#' @family accessors
#' @export
contrast <- function(fit, method, ci) {
  stopifnot(inherits(fit, "separable_effects"))

  if (missing(ci) || is.null(ci)) {
    stop(
      "'ci' argument is required. Contrasts without uncertainty are not ",
      "reported.\n  Compute bootstrap first: boot <- bootstrap(fit, n_boot = 500)\n",
      "  Then: contrast(fit, method = '<name>', ci = boot)",
      call. = FALSE
    )
  }
  stopifnot(inherits(ci, "causal_cr_bootstrap"))

  if (missing(method)) {
    stop(
      "'method' argument is required. Available: ",
      paste(names(fit$cumulative_incidence), collapse = ", "),
      call. = FALSE
    )
  }

  contrasts_df <- compute_contrasts(fit, method = method, ci = ci)

  structure(
    list(
      contrasts = contrasts_df,
      method    = method,
      alpha     = ci$alpha
    ),
    class = "causal_cr_contrast"
  )
}


#' Extract Diagnostics
#'
#' Returns a `"causal_cr_diagnostic"` object combining weight-level
#' diagnostics (from IPW, when applicable) with model-level diagnostics
#' (from all fitted hazard models). Has its own `print()` and `plot()`
#' methods.
#'
#' @param fit A `"separable_effects"` object from [separable_effects()].
#'
#' @return An S3 object of class `"causal_cr_diagnostic"` with:
#'   \describe{
#'     \item{weight_summary}{Data frame summarizing weight distributions
#'       (raw and truncated). NULL if no IPW method was run.}
#'     \item{flagged_ids}{IDs affected by weight truncation. Empty if no IPW.}
#'     \item{truncate}{Length-2 percentile bounds used, or NULL.}
#'     \item{model_checks}{Named list (`y`, `d`, `c`) of per-model
#'       diagnostics: convergence, min/max fitted probabilities, positivity
#'       violation flag, captured glm warnings.}
#'   }
#'
#' @seealso [risk()], [contrast()], [plot.causal_cr_diagnostic()]
#' @family accessors
#' @export
diagnostic <- function(fit) {
  stopifnot(inherits(fit, "separable_effects"))

  weights_slot <- fit$weights

  out <- list(
    weight_summary = weights_slot$weight_summary,
    flagged_ids    = weights_slot$flagged_ids %||% integer(0),
    flagged_log    = weights_slot$flagged_log,
    truncate       = weights_slot$truncate,
    model_checks   = fit$model_checks
  )

  structure(out, class = "causal_cr_diagnostic")
}
