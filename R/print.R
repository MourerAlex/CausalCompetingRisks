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
#' Prints per-method cumulative incidence at the final time, a brief model
#' checks summary, and guidance toward the proper accessors for contrasts
#' and confidence intervals.
#'
#' @param object A `"separable_effects"` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the per-method list of cumulative incidence
#'   data frames.
#' @export
summary.separable_effects <- function(object, ...) {
  cat("Separable Effects - separable_effects\n")
  cat("=============================\n\n")
  cat("Method(s):", paste(names(object$cumulative_incidence),
                          collapse = ", "),
      "| N:", object$n, "\n")
  cat("Time points:", length(object$times), "\n\n")

  cat("Cumulative incidence at final time:\n")
  for (m in names(object$cumulative_incidence)) {
    df <- object$cumulative_incidence[[m]]
    last <- df[nrow(df), ]
    cat(sprintf(
      "  [%s]  (1,1)=%.4f  (0,0)=%.4f  (1,0)=%.4f  (0,1)=%.4f\n",
      m, last$arm_11, last$arm_00, last$arm_10, last$arm_01
    ))
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

  cat("\nFor causal contrasts with confidence intervals:\n")
  cat("  boot <- bootstrap(fit, n_boot = 500)\n")
  cat("  contrast(fit, method = '<name>', ci = boot)\n")

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


#' Print a causal_cr_contrast Object
#'
#' Shows the contrast table at the final time point (compact view), plus
#' the method, significance level, and a pointer to the full long-format
#' data frame.
#'
#' @param x A `"causal_cr_contrast"` object from [contrast()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_cr_contrast <- function(x, ...) {
  cat("Causal Contrasts (causal_cr_contrast)\n")
  cat("--------------------------------------\n")
  cat("Method: ", x$method, "\n", sep = "")
  cat(sprintf("Significance level: %g (%.0f%% CIs)\n\n",
              x$alpha, (1 - x$alpha) * 100))

  # Compact view: contrasts at the final time point only
  final_k <- max(x$contrasts$k)
  at_final <- x$contrasts[x$contrasts$k == final_k, ]

  cat(sprintf("At final time (k = %g):\n", final_k))
  print(
    at_final[, c("contrast", "decomp", "measure", "estimate", "lower", "upper")],
    row.names = FALSE
  )

  cat(sprintf(
    "\nFull long-format data: %d rows (%d time points x 10 contrasts). Use x$contrasts.\n",
    nrow(x$contrasts),
    length(unique(x$contrasts$k))
  ))
  invisible(x)
}


#' Print a causal_cr_diagnostic Object
#'
#' Shows per-model fit diagnostics (convergence, fitted probability range,
#' positivity violation flag) and IPW weight summary if available. Uses
#' the `flagged_ids` / `flagged_log` fields (both truncation and trimming
#' cases).
#'
#' @param x A `"causal_cr_diagnostic"` object from [diagnostic()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_cr_diagnostic <- function(x, ...) {
  cat("Diagnostics (causal_cr_diagnostic)\n")
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
