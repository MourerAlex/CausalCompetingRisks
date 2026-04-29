#' Extract Cumulative Incidence Curves
#'
#' Returns a `"separable_effects_risk"` object containing cumulative
#' incidence estimates for each treatment arm in long format. Optionally
#' pairs with a bootstrap object to fold confidence bands into the same
#' table.
#'
#' @param fit A `"separable_effects"` object from [separable_effects()].
#' @param ci Optional. A `"separable_effects_bootstrap"` object from
#'   [bootstrap()]. When provided, the `lower` / `upper` columns of
#'   `$risk` are populated; otherwise they are `NA_real_`.
#'
#' @return An S3 object of class `"separable_effects_risk"` with:
#'   \describe{
#'     \item{risk}{Long-format data.frame, one row per
#'       `(method, arm, k)` triple, with columns `method`, `arm`,
#'       `a_y`, `a_d`, `k`, `value`, `lower`, `upper`. Sorted by
#'       method, arm (per [arm_spec()] order), then time.}
#'     \item{replicates}{4D bootstrap array (or NULL). Kept for
#'       per-contrast CIs in [plot.separable_effects_risk()].}
#'     \item{alpha}{Bootstrap significance level (or NULL).}
#'     \item{person_time, id_col, treatment_col, times}{References to
#'       the underlying fit (used by the optional risk-table panel
#'       inside the plot method).}
#'   }
#'
#' @seealso [contrast()], [diagnostic()], [plot.separable_effects_risk()],
#'   [bootstrap()]
#' @family accessors
#' @export
risk <- function(fit, ci = NULL) {
  stopifnot(inherits(fit, "separable_effects"))
  if (!is.null(ci)) {
    stopifnot(inherits(ci, "separable_effects_bootstrap"))
  }
  structure(
    list(
      risk          = build_risk_long(fit, ci),
      # Full bootstrap replicates + alpha for PROPER per-contrast CIs in
      # plot.separable_effects_risk (per-replicate difference, then quantile).
      replicates    = if (!is.null(ci)) ci$replicates else NULL,
      alpha         = if (!is.null(ci)) ci$alpha      else NULL,
      # References needed for plot's risk_table option (avoids forcing
      # the user to pass `fit` again at plot time). These are R pointers,
      # not copies — no real memory overhead unless the underlying data
      # is modified.
      person_time   = fit$person_time,
      id_col        = fit$id_col,
      treatment_col = fit$treatment_col,
      times         = fit$times
    ),
    class = "separable_effects_risk"
  )
}


#' Build the Long-Format `$risk` Data.frame
#'
#' Pivots `fit$cumulative_incidence` (per-method wide) and (optionally)
#' `ci$ci_curves` (per-method wide bands) into a single long table
#' `(method, arm, a_y, a_d, k, value, lower, upper)`. The `(a_y, a_d)`
#' columns come from [arm_spec()].
#'
#' @param fit A `"separable_effects"` object.
#' @param ci A `"separable_effects_bootstrap"` object or NULL.
#' @return Long-format data.frame; `lower` / `upper` are `NA_real_` when
#'   `ci` is NULL or for arms missing from a method's bands.
#' @family internal
#' @keywords internal
build_risk_long <- function(fit, ci) {
  spec <- arm_spec()
  cum_inc_list <- fit$cumulative_incidence
  ci_curves    <- if (!is.null(ci)) ci$ci_curves else NULL

  rows <- list()
  for (m in names(cum_inc_list)) {
    wide <- cum_inc_list[[m]]
    bands <- if (!is.null(ci_curves)) ci_curves[[m]] else NULL
    for (i in seq_len(nrow(spec))) {
      arm <- spec$name[i]
      if (!arm %in% names(wide)) next  # method emitted a strict subset
      lo <- if (!is.null(bands) && paste0(arm, "_lower") %in% names(bands)) {
        bands[[paste0(arm, "_lower")]]
      } else {
        rep(NA_real_, length(wide$k))
      }
      hi <- if (!is.null(bands) && paste0(arm, "_upper") %in% names(bands)) {
        bands[[paste0(arm, "_upper")]]
      } else {
        rep(NA_real_, length(wide$k))
      }
      rows[[length(rows) + 1L]] <- data.frame(
        method = m,
        arm    = arm,
        a_y    = spec$a_y[i],
        a_d    = spec$a_d[i],
        k      = wide$k,
        value  = wide[[arm]],
        lower  = lo,
        upper  = hi,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      method = character(), arm = character(),
      a_y = integer(), a_d = integer(),
      k = numeric(), value = numeric(),
      lower = numeric(), upper = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}


#' Extract Causal Contrasts (Long Format with Confidence Intervals)
#'
#' Returns a `"separable_effects_contrast"` object with a long-format data frame
#' of total, separable direct, and separable indirect effects for a single
#' estimation method. Bootstrap confidence intervals are required — the
#' package does not report contrasts without uncertainty.
#'
#' @param fit A `"separable_effects"` object from [separable_effects()].
#' @param method Character (length 1). Which method's cumulative incidence
#'   to contrast. Must be in `names(fit$cumulative_incidence)`.
#' @param ci A `"separable_effects_bootstrap"` object from [bootstrap()]. Required.
#'
#' @return An S3 object of class `"separable_effects_contrast"` with:
#'   \describe{
#'     \item{contrasts}{Long-format data frame with columns `k`, `contrast`,
#'       `decomp`, `measure`, `estimate`, `lower`, `upper`. 10 rows per
#'       time point (1 total + 2 direct + 2 indirect, RD and RR each).}
#'     \item{method}{The method used.}
#'     \item{alpha}{Significance level from the bootstrap object.}
#'   }
#'
#' @seealso [risk()], [diagnostic()], [plot.separable_effects_contrast()],
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
  stopifnot(inherits(ci, "separable_effects_bootstrap"))

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
    class = "separable_effects_contrast"
  )
}


#' Extract Diagnostics
#'
#' Returns a `"separable_effects_diagnostic"` object combining weight-level
#' diagnostics (from IPW, when applicable) with model-level diagnostics
#' (from all fitted hazard models). Has its own `print()` and `plot()`
#' methods.
#'
#' @param fit A `"separable_effects"` object from [separable_effects()].
#'
#' @return An S3 object of class `"separable_effects_diagnostic"` with:
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
#' @seealso [risk()], [contrast()], [plot.separable_effects_diagnostic()]
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

  structure(out, class = "separable_effects_diagnostic")
}
