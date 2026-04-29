#' Print a separable_effects Object
#'
#' @param x A `"separable_effects"` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns `x`.
#' @export
print.separable_effects <- function(x, ...) {
  cat("Separable Effects Estimation\n")
  cat("----------------------------\n")

  # One-line identification header. Surfaces the estimand framing only;
  # the swap-weight numeric read-out is intentionally suppressed pending
  # a math check against the Stensrud Appendix — see `@section Isolation
  # read-out (deferred)` on [assumptions()] for context. To restore, see
  # commit history.
  cat("Estimand: P(Y^{a_Y, a_D, c_bar=0}_{K+1}=1)",
      "under separable-effects identification\n")
  cat("See `assumptions(fit)` for the full identification block.\n")

  cat("Method(s):", paste(names(x$cumulative_incidence), collapse = ", "), "\n")
  cat("N subjects:", x$n, "\n")
  cat("Time points:", length(x$times), "\n")

  # Per-method cumulative incidence at final time
  cat("\nCumulative incidence at final time:\n")
  for (m in names(x$cumulative_incidence)) {
    df <- x$cumulative_incidence[[m]]
    last <- df[nrow(df), ]
    cat(sprintf(
      "  [%s]  (1,1)=%.4f  (0,0)=%.4f  (1,0)=%.4f  (0,1)=%.4f\n",
      m, last$arm_11, last$arm_00, last$arm_10, last$arm_01
    ))
  }

  # Model-checks summary (no binary positivity flag — just count the
  # continuous signals we surface)
  if (!is.null(x$model_checks)) {
    issues <- 0
    for (chk in x$model_checks) {
      if (is.null(chk)) next
      if (!isTRUE(chk$converged)) issues <- issues + 1
      if (length(chk$glm_warnings) > 0) issues <- issues + 1
    }
    if (issues > 0) {
      cat("\nModel checks:", issues,
          "issue(s) - use `fit$model_checks` to inspect.\n")
    }
  }

  if (length(x$warnings) > 0) {
    cat("\nWarnings:", length(x$warnings),
        "- use `fit$warnings` to inspect.\n")
  }

  cat("\nUse risk(), contrast(), diagnostic() to extract components.\n")
  invisible(x)
}


#' Summary of a separable_effects Object
#'
#' Prints per-method cumulative incidence at the selected time, a brief
#' model checks summary, and (when a bootstrap is supplied) the contrast
#' table at that same time.
#'
#' @param object A `"separable_effects"` object.
#' @param ci Optional. A `"separable_effects_bootstrap"` object from
#'   [bootstrap()]. When provided, the contrast table at the selected
#'   `time` is printed alongside the cumulative incidence summary.
#' @param time Numeric scalar or NULL. Time point at which to summarise.
#'   NULL (default) selects the final cut time (`max(object$times)`). A
#'   user-supplied value is snapped to the nearest cut time and a
#'   `message()` is emitted when snapping changes the value.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the per-method list of cumulative incidence
#'   data frames.
#' @export
summary.separable_effects <- function(object, ci = NULL, time = NULL, ...) {
  if (!is.null(ci)) {
    stopifnot(inherits(ci, "separable_effects_bootstrap"))
  }

  k_at <- snap_time(time, object$times)

  cat("Separable Effects - separable_effects\n")
  cat("=============================\n\n")
  cat("Method(s):", paste(names(object$cumulative_incidence),
                          collapse = ", "),
      "| N:", object$n, "\n")
  cat("Time points:", length(object$times), "\n\n")

  cat(sprintf("Cumulative incidence at k = %g:\n", k_at))
  for (m in names(object$cumulative_incidence)) {
    df <- object$cumulative_incidence[[m]]
    row <- df[df$k == k_at, , drop = FALSE]
    if (nrow(row) == 0L) next
    cat(sprintf(
      "  [%s]  (1,1)=%.4f  (0,0)=%.4f  (1,0)=%.4f  (0,1)=%.4f\n",
      m, row$arm_11, row$arm_00, row$arm_10, row$arm_01
    ))
  }

  # Contrast block (only when bootstrap is attached)
  if (!is.null(ci)) {
    cat(sprintf(
      "\nContrasts at k = %g (%.0f%% CIs):\n", k_at, (1 - ci$alpha) * 100
    ))
    for (m in names(object$cumulative_incidence)) {
      if (!m %in% dimnames(ci$replicates)[[2]]) next
      ctr_full <- compute_contrasts(object, method = m, ci = ci)
      ctr_at   <- ctr_full[ctr_full$k == k_at &
                             ctr_full$measure == "rd", , drop = FALSE]
      cat(sprintf("  [%s]\n", m))
      for (i in seq_len(nrow(ctr_at))) {
        r <- ctr_at[i, ]
        decomp_label <- if (is.na(r$decomp)) ""
                        else paste0(" (", r$decomp, ")")
        cat(sprintf(
          "    %-8s%-4s  RD = %6.3f  [%6.3f, %6.3f]\n",
          r$contrast, decomp_label,
          r$estimate, r$lower, r$upper
        ))
      }
    }
  }

  # Model-checks line — report non-convergence. Positivity is a continuous
  # signal; users should inspect diagnostic(fit) for min_fitted / max_fitted.
  if (!is.null(object$model_checks)) {
    non_converged <- 0L
    for (chk in object$model_checks) {
      if (is.null(chk)) next
      if (!isTRUE(chk$converged)) non_converged <- non_converged + 1L
    }
    if (non_converged > 0) {
      cat(sprintf(
        "\nModel checks: %d non-converged model(s). See `fit$model_checks`.\n",
        non_converged
      ))
    }
  }

  if (is.null(ci)) {
    cat("\nFor causal contrasts with confidence intervals:\n")
    cat("  boot <- bootstrap(fit, n_boot = 500)\n")
    cat("  summary(fit, ci = boot)        # contrasts at final time\n")
    cat("  contrast(fit, method = '<name>', ci = boot)\n")
  }

  invisible(object$cumulative_incidence)
}


#' Confidence Intervals for separable_effects
#'
#' Deprecated pathway. In the current design, confidence intervals are
#' computed by pairing a fit with a separate [bootstrap()] object.
#'
#' @param object A `"separable_effects"` object.
#' @param parm Not used (included for S3 consistency).
#' @param level Not used.
#' @param ... Additional arguments (currently unused).
#'
#' @return This method always errors — use the `bootstrap()` + `contrast()`
#'   pattern instead.
#' @export
confint.separable_effects <- function(object, parm = NULL, level = 0.95, ...) {
  stop(
    "Confidence intervals are not stored inside `fit` in this package. ",
    "Compute them explicitly:\n",
    "  boot <- bootstrap(fit, n_boot = 500)\n",
    "  contrast(fit, method = '<name>', ci = boot)",
    call. = FALSE
  )
}


#' Print a separable_effects_risk Object
#'
#' Shows per-method cumulative incidence at the final time point, indicates
#' whether bootstrap CI bands are attached, and points toward `plot()`.
#'
#' @param x A `"separable_effects_risk"` object from [risk()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.separable_effects_risk <- function(x, ...) {
  cat("Cumulative Incidence Curves (separable_effects_risk)\n")
  cat("---------------------------------------------\n")

  methods_avail <- names(x$cumulative_incidence)
  cat("Methods: ", paste(methods_avail, collapse = ", "), "\n", sep = "")
  cat("Bootstrap CIs: ", if (!is.null(x$ci_curves)) "yes" else "no", "\n\n",
      sep = "")

  for (m in methods_avail) {
    df <- x$cumulative_incidence[[m]]
    last <- df[nrow(df), ]
    cat(sprintf(
      "[%s] at final time (k = %g):\n  (1,1)=%.4f  (0,0)=%.4f  (1,0)=%.4f  (0,1)=%.4f\n\n",
      m, df$k[nrow(df)],
      last$arm_11, last$arm_00, last$arm_10, last$arm_01
    ))
  }

  cat("Use plot(risk(fit), method = '<name>') to visualize.\n")
  invisible(x)
}


#' Print a separable_effects_contrast Object
#'
#' Shows the contrast table at the selected time point, plus the method
#' and significance level.
#'
#' @param x A `"separable_effects_contrast"` object from [contrast()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.separable_effects_contrast <- function(x, ...) {
  cat("Causal Contrasts (separable_effects_contrast)\n")
  cat("--------------------------------------\n")
  cat("Method: ", x$method, "\n", sep = "")
  cat(sprintf("Significance level: %g (%.0f%% CIs)\n\n",
              x$alpha, (1 - x$alpha) * 100))

  k_at <- x$time %||% max(x$contrasts$k)
  cat(sprintf("At k = %g:\n", k_at))
  print(
    x$contrasts[, c("contrast", "decomp", "measure", "estimate", "lower", "upper")],
    row.names = FALSE
  )

  cat(sprintf(
    "\n%d rows (1 time point x 10 contrasts). Pass `time = ...` to ",
    nrow(x$contrasts)
  ))
  cat("contrast() for a different cut time.\n")
  invisible(x)
}


#' Print a separable_effects_diagnostic Object
#'
#' Shows per-model fit diagnostics (convergence, fitted probability range,
#' positivity violation flag) and IPW weight summary if available. Uses
#' the `flagged_ids` / `flagged_log` fields (both truncation and trimming
#' cases).
#'
#' @param x A `"separable_effects_diagnostic"` object from [diagnostic()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.separable_effects_diagnostic <- function(x, ...) {
  cat("Diagnostics (separable_effects_diagnostic)\n")
  cat("-----------------------------------\n")

  # --- Model checks ---
  if (!is.null(x$model_checks)) {
    cat("\nModel checks:\n")
    any_printed <- FALSE
    for (m in names(x$model_checks)) {
      chk <- x$model_checks[[m]]
      if (is.null(chk)) next
      any_printed <- TRUE
      cat(sprintf(
        "  [%s]  converged=%s  fitted in [%.3g, %.3g]\n",
        chk$label %||% m,
        if (isTRUE(chk$converged)) "TRUE" else "FALSE",
        chk$min_fitted, chk$max_fitted
      ))
      if (length(chk$glm_warnings) > 0) {
        cat("    glm warnings:\n")
        for (w in chk$glm_warnings) {
          cat("      - ", w, "\n", sep = "")
        }
      }
    }
    if (!any_printed) cat("  (no fitted models)\n")
  }

  # --- Weight summary (IPW only) ---
  if (!is.null(x$weight_summary)) {
    cat("\nWeight summary:\n")
    print(x$weight_summary, row.names = FALSE)

    n_flagged <- length(x$flagged_ids)
    n_log_rows <- if (!is.null(x$flagged_log)) nrow(x$flagged_log) else 0
    trunc_label <- if (is.null(x$truncate)) "none" else
      sprintf("[%g, %g]", x$truncate[1], x$truncate[2])
    cat(sprintf(
      "\nFlagged subjects (truncate=%s): %d  (%d row(s) in flagged_log)\n",
      trunc_label, n_flagged, n_log_rows
    ))
  } else {
    cat("\nNo IPW weight diagnostics (no IPW method ran).\n")
  }

  invisible(x)
}
