#' Print a causal_cr Object
#'
#' @param x A `"causal_cr"` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns `x`.
#' @export
print.causal_cr <- function(x, ...) {
  cat("Separable Effects Estimation\n")
  cat("----------------------------\n")
  cat("Method:", x$method, "\n")
  cat("N subjects:", x$n, "\n")
  cat("Time points:", length(x$times), "\n")
  cat("Bootstrap:", if (!is.null(x$bootstrap)) "Yes" else "No", "\n")
  cat("\nCumulative incidence at final time:\n")
  last <- x$cumulative_incidence[nrow(x$cumulative_incidence), ]
  cat("  arm(1,1):", round(last$arm_11, 4), "\n")
  cat("  arm(0,0):", round(last$arm_00, 4), "\n")
  cat("  arm(1,0):", round(last$arm_10, 4), "\n")
  cat("\nUse risk(), contrast(), diagnostic() to extract components.\n")
  invisible(x)
}


#' Summary of a causal_cr Object
#'
#' @param object A `"causal_cr"` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the contrasts data frame.
#' @export
summary.causal_cr <- function(object, ...) {
  cat("Separable Effects — Causal Contrasts\n")
  cat("=====================================\n\n")
  cat("Method:", object$method, "| N:", object$n, "\n\n")

  ct <- object$contrasts
  last <- ct[nrow(ct), ]

  cat("At final time (k =", last$k, "):\n\n")
  cat(sprintf("  %-25s  RD = %+.4f   RR = %.4f\n",
              "Total effect", last$total_rd, last$total_rr))
  cat(sprintf("  %-25s  RD = %+.4f   RR = %.4f\n",
              "Sep. direct (A_Y)", last$sep_direct_rd, last$sep_direct_rr))
  cat(sprintf("  %-25s  RD = %+.4f   RR = %.4f\n",
              "Sep. indirect (A_D)", last$sep_indirect_rd, last$sep_indirect_rr))

  if (!is.null(object$bootstrap)) {
    cat("\nBootstrap 95% CI available. Use confint() for details.\n")
  }

  invisible(ct)
}


#' Confidence Intervals for causal_cr
#'
#' Prints bootstrap confidence intervals for causal contrasts at eval_times.
#'
#' @param object A `"causal_cr"` object with bootstrap results.
#' @param parm Not used (included for S3 consistency).
#' @param level Not used (alpha is set at estimation time).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the CI data frame.
#' @export
confint.causal_cr <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(object$bootstrap)) {
    stop("No bootstrap results. Re-run causal_cr() with bootstrap = TRUE.",
         call. = FALSE)
  }

  cat("Bootstrap Confidence Intervals\n")
  cat("==============================\n\n")
  print(object$bootstrap$ci_contrasts)

  invisible(object$bootstrap$ci_contrasts)
}


#' Print a causal_cr_risk Object
#'
#' @param x A `"causal_cr_risk"` object from [risk()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_cr_risk <- function(x, ...) {
  cat("Cumulative Incidence Curves\n")
  cat("---------------------------\n")
  print(utils::head(x$cumulative_incidence, 10))
  if (nrow(x$cumulative_incidence) > 10) {
    cat("... (", nrow(x$cumulative_incidence), " time points total)\n")
  }
  invisible(x)
}


#' Print a causal_cr_contrast Object
#'
#' @param x A `"causal_cr_contrast"` object from [contrast()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_cr_contrast <- function(x, ...) {
  cat("Causal Contrasts (Separable Effects)\n")
  cat("-------------------------------------\n")
  print(utils::head(x$contrasts, 10))
  if (nrow(x$contrasts) > 10) {
    cat("... (", nrow(x$contrasts), " time points total)\n")
  }
  invisible(x)
}


#' Print a causal_cr_diagnostic Object
#'
#' @param x A `"causal_cr_diagnostic"` object from [diagnostic()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_cr_diagnostic <- function(x, ...) {
  cat("IPW Diagnostics\n")
  cat("----------------\n")
  if (!is.null(x$weight_summary)) {
    cat("\nWeight summary:\n")
    print(x$weight_summary)
    cat("\nTruncated/trimmed IDs:", length(x$truncated_ids), "\n")
  } else {
    cat("No IPW diagnostics available.\n")
  }
  invisible(x)
}
