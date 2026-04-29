#' Identifying Assumptions for a Separable-Effects Fit
#'
#' Returns a structured object listing the target estimand and the
#' identifying assumptions under which it is identified from the observed
#' data, plus a numeric isolation read-out (no verdict). Designed for
#' inspection at the console (`print()`) and reuse in vignettes / papers
#' (`format(x, style = "markdown")`).
#'
#' The accessor does not run any new computation on the fit beyond
#' summarising the per-row swap weights for the isolation slot. All
#' statements are textual and refer to the estimand and the conditions
#' formalised in Stensrud et al. (2020, 2021).
#'
#' @param fit A `"separable_effects"` object from [separable_effects()].
#'
#' @return An S3 object of class `"separable_effects_assumptions"` with:
#'   \describe{
#'     \item{estimand}{Named list. Target counterfactual notation, the
#'       four arms `(a_Y, a_D)`, and the two decompositions (A and B)
#'       expressed as differences between arm cumulative incidence.}
#'     \item{assumptions}{Named list with one entry per identifying
#'       assumption (GDA, treatment exchangeability, censoring
#'       exchangeability, consistency, positivity). Each entry is a
#'       list with fields `statement`, `testable`, `citation`,
#'       `diagnostic`.}
#'     \item{dismissible_D1}{Named list (statement / testable / citation
#'       / diagnostic) for the dismissible component condition Δ1
#'       (cause-specific hazard of Y does not depend on `a_D`).}
#'     \item{dismissible_D2}{Named list for Δ2 (cause-specific hazard
#'       of D does not depend on `a_Y`).}
#'   }
#'
#' @section Isolation read-out (deferred):
#' The locked API spec calls for a numeric isolation read-out — per-row
#' swap-weight ranges `w_d` and `w_y`, no verdict — both in the
#' `$isolation` slot of this object and as a one-line header on
#' `print(fit)`. The empirical mapping `w_d ≈ 1 ↔ Δ2` (and `w_y ≈ 1 ↔
#' Δ1`) has not yet been math-checked against the Stensrud Appendix.
#' Until that check is done, the read-out is intentionally omitted to
#' avoid surfacing an unverified causal claim. The textual statements
#' for Δ1 and Δ2 (with diagnostic pointers) remain in place. To
#' restore: re-include the `$isolation` slot built by
#' [isolation_summary()] (kept for that purpose) and the corresponding
#' sections in `print()` / `format()` and the one-line header in
#' [print.separable_effects()]; see commit history for the prior
#' implementation.
#'
#' @references
#' Stensrud MJ, Young JG, Didelez V, Robins JM, Hernán MA (2020).
#' Separable Effects for Causal Inference in the Presence of Competing
#' Events. \doi{10.1080/01621459.2020.1765783}
#'
#' Stensrud MJ, Hernán MA, Tchetgen Tchetgen EJ, Robins JM, Didelez V,
#' Young JG (2021). A generalized theory of separable effects in
#' competing event settings. \doi{10.1007/s10985-021-09530-8}
#'
#' @seealso [separable_effects()], [risk()], [contrast()], [diagnostic()]
#' @family accessors
#' @export
assumptions <- function(fit) {
  stopifnot(inherits(fit, "separable_effects"))

  estimand <- list(
    target = "P(Y^{a_Y, a_D, c_bar = 0}_{K+1} = 1)",
    description = paste0(
      "Cumulative incidence of the primary event Y by interval K+1 ",
      "under hypothetical assignment to the four-arm regime ",
      "(a_Y, a_D), with censoring eliminated."
    ),
    arms = c("(1,1)", "(0,0)", "(1,0)", "(0,1)"),
    decomposition_A = list(
      label = "Decomposition A (vary A_Y first; intermediate arm = arm_10)",
      sde = "SDE-A = F^(1,0) - F^(0,0) = effect of A_Y at a_D = 0",
      sie = "SIE-A = F^(1,1) - F^(1,0) = effect of A_D at a_Y = 1"
    ),
    decomposition_B = list(
      label = "Decomposition B (vary A_D first; intermediate arm = arm_01)",
      sie = "SIE-B = F^(0,1) - F^(0,0) = effect of A_D at a_Y = 0",
      sde = "SDE-B = F^(1,1) - F^(0,1) = effect of A_Y at a_D = 1"
    )
  )

  assumptions_list <- list(

    GDA = list(
      statement = paste0(
        "Generalized decomposition: the observed treatment A can be ",
        "decomposed into components A_Y and A_D such that A = A_Y = A_D ",
        "in the observed data. In a future study, A_Y and A_D could in ",
        "principle be assigned different values; setting both to the ",
        "same value recovers the observed treatment."
      ),
      testable = paste0(
        "Untestable in a two-arm trial. Testable in a hypothetical ",
        "four-arm trial that separately randomises A_Y and A_D."
      ),
      citation = "Stensrud et al. (2020), Section 3; Stensrud et al. (2021).",
      diagnostic = NULL
    ),

    E1_treatment_exchangeability = list(
      statement = paste0(
        "Counterfactual outcomes under no censoring are independent of ",
        "the assigned treatment given measured covariates: ",
        "(Y^{a, c_bar=0}, D^{a, c_bar=0}) _||_ A | L."
      ),
      testable = paste0(
        "Holds by design in a randomised controlled trial. In ",
        "observational data, untestable from the observed marginals; ",
        "covariate balance in propensity-score diagnostics is supportive ",
        "but not conclusive."
      ),
      citation = "Stensrud et al. (2020), Section 3.2.",
      diagnostic = NULL
    ),

    E2_censoring_exchangeability = list(
      statement = paste0(
        "Counterfactual outcomes under no censoring are independent of ",
        "the censoring event given history: ",
        "(Y^{a, c_bar=0}_{k+1}, D^{a, c_bar=0}_{k+1}) _||_ C_{k+1} | ",
        "Y_k = D_k = C_bar_k = 0, L, A. Equivalent to ignorable censoring."
      ),
      testable = "Untestable from observed data alone.",
      citation = "Stensrud et al. (2020), Section 3.3.",
      diagnostic = paste0(
        "When `censoring_weights = TRUE`, the IPCW estimator addresses ",
        "this assumption. The censoring model fit is summarised in ",
        "diagnostic(fit)$model_checks$c."
      )
    ),

    consistency = list(
      statement = paste0(
        "If A = a and C_bar_k = 0, then the observed event indicators ",
        "equal their counterfactual analogues: Y^{a, c_bar=0}_k = Y_k ",
        "and D^{a, c_bar=0}_k = D_k. Applies only in the observed arms ",
        "(a_Y = a_D) and only on uncensored rows."
      ),
      testable = "Untestable; assumed by design when treatment is well-defined.",
      citation = "Stensrud et al. (2020), Section 3.4.",
      diagnostic = NULL
    ),

    positivity = list(
      statement = paste0(
        "For every covariate pattern, both treatments occur with ",
        "positive probability, and uncensored observation persists at ",
        "every interval where the joint of (A, history, L) has positive ",
        "support. Without positivity, the identifying formulas divide ",
        "by zero."
      ),
      testable = paste0(
        "Diagnosable in fitted-probability ranges of the hazard and ",
        "propensity models. Probabilities very close to 0 or 1 indicate ",
        "near-violations."
      ),
      citation = "Stensrud et al. (2020), Section 3.5.",
      diagnostic = paste0(
        "diagnostic(fit)$model_checks$<m>$min_fitted / $max_fitted for ",
        "each fitted model m in {y, d, c, a}."
      )
    )
  )

  dismissible_D1 <- list(
    statement = paste0(
      "Δ1 — dismissible component condition for Y. The cause-specific ",
      "hazard of Y among event-free individuals does not depend on ",
      "a_D, conditional on the covariate history."
    ),
    testable = paste0(
      "Untestable from observed data. A biological / causal claim, not ",
      "a statistical property of the data."
    ),
    citation = paste0(
      "Stensrud et al. (2020), Section 3.6 (Δ1); Stensrud et al. ",
      "(2021), equation (34) for the time-varying form."
    ),
    diagnostic = paste0(
      "Empirical numeric read-out (swap-weight ranges for w_y) is ",
      "deferred pending a math check against the Stensrud Appendix; ",
      "see `?assumptions` for context."
    )
  )

  dismissible_D2 <- list(
    statement = paste0(
      "Δ2 — dismissible component condition for D. The cause-specific ",
      "hazard of D among event-free individuals does not depend on ",
      "a_Y, conditional on the covariate history."
    ),
    testable = "Untestable from observed data.",
    citation = paste0(
      "Stensrud et al. (2020), Section 3.6 (Δ2); Stensrud et al. ",
      "(2021), equation (35) for the time-varying form."
    ),
    diagnostic = paste0(
      "Empirical numeric read-out (swap-weight ranges for w_d) is ",
      "deferred pending a math check against the Stensrud Appendix; ",
      "see `?assumptions` for context."
    )
  )

  # --- Isolation slot intentionally omitted; see @section above. ---
  # The helper isolation_summary() is kept for the future re-enable.

  structure(
    list(
      estimand        = estimand,
      assumptions     = assumptions_list,
      dismissible_D1  = dismissible_D1,
      dismissible_D2  = dismissible_D2
    ),
    class = "separable_effects_assumptions"
  )
}


#' Compute Isolation Summary (No Verdict) — currently unused
#'
#' Pulls the per-row swap weights from the fit and reports their numeric
#' range. Returns NULL fields if no IPW method ran (i.e., gformula-only).
#'
#' Currently NOT consumed by [assumptions()] or [print.separable_effects()];
#' kept for the future re-enable once the empirical link
#' `w_d ≈ 1 ↔ Δ2` (and `w_y ≈ 1 ↔ Δ1`) is math-checked against the
#' Stensrud Appendix. See `@section Isolation read-out (deferred)` on
#' [assumptions()].
#'
#' @param fit A `"separable_effects"` object.
#' @return A named list with `range_w_d`, `range_w_y`, `n_rows_w_d`,
#'   `n_rows_w_y`, and `note`.
#' @family internal
#' @keywords internal
isolation_summary <- function(fit) {

  pt <- fit$weights$pt_data_weighted

  range_or_null <- function(cols) {
    if (is.null(pt)) return(list(range = NULL, n = 0L))
    present <- intersect(cols, names(pt))
    if (length(present) == 0L) return(list(range = NULL, n = 0L))
    vals <- unlist(lapply(present, function(co) pt[[co]]), use.names = FALSE)
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) return(list(range = NULL, n = 0L))
    list(range = range(vals), n = length(vals))
  }

  rd <- range_or_null(c("w_d_arm_10", "w_d_arm_01"))
  ry <- range_or_null(c("w_y_arm_10", "w_y_arm_01"))

  note <- paste0(
    "Isolation cannot be falsified from data. Component hazard ratios ",
    "(swap weights) near 1 are consistent with — but do not prove — Δ1 ",
    "(for w_y) and Δ2 (for w_d). Sensitivity analysis under partial ",
    "isolation is the appropriate response when the ranges depart ",
    "meaningfully from 1; see Stensrud & Young (2021)."
  )

  list(
    range_w_d  = rd$range,
    range_w_y  = ry$range,
    n_rows_w_d = rd$n,
    n_rows_w_y = ry$n,
    note       = note
  )
}


#' Print a separable_effects_assumptions Object
#'
#' Renders a human-readable identification block: estimand, the standard
#' identifying assumptions, and the two dismissible component conditions
#' Δ1 and Δ2. Numeric isolation read-out is currently deferred; see
#' `@section Isolation read-out (deferred)` on [assumptions()].
#'
#' @param x A `"separable_effects_assumptions"` object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.separable_effects_assumptions <- function(x, ...) {

  cat("Identifying Assumptions (separable_effects)\n")
  cat("-------------------------------------------\n\n")

  # --- Estimand ---
  cat("Estimand:\n")
  cat("  ", x$estimand$target, "\n", sep = "")
  cat("  ", x$estimand$description, "\n", sep = "")
  cat("\nArms (a_Y, a_D): ",
      paste(x$estimand$arms, collapse = ", "), "\n", sep = "")
  cat("\n", x$estimand$decomposition_A$label, "\n", sep = "")
  cat("  ", x$estimand$decomposition_A$sde, "\n", sep = "")
  cat("  ", x$estimand$decomposition_A$sie, "\n", sep = "")
  cat("\n", x$estimand$decomposition_B$label, "\n", sep = "")
  cat("  ", x$estimand$decomposition_B$sie, "\n", sep = "")
  cat("  ", x$estimand$decomposition_B$sde, "\n", sep = "")

  # --- Assumptions ---
  cat("\nIdentifying assumptions:\n")
  for (nm in names(x$assumptions)) {
    a <- x$assumptions[[nm]]
    cat("\n  [", nm, "]\n", sep = "")
    cat("    ", a$statement, "\n", sep = "")
    cat("    Testable: ", a$testable, "\n", sep = "")
    cat("    Citation: ", a$citation, "\n", sep = "")
    if (!is.null(a$diagnostic)) {
      cat("    Diagnostic: ", a$diagnostic, "\n", sep = "")
    }
  }

  # --- Dismissible component conditions ---
  for (nm in c("dismissible_D1", "dismissible_D2")) {
    a <- x[[nm]]
    cat("\n  [", nm, "]\n", sep = "")
    cat("    ", a$statement, "\n", sep = "")
    cat("    Testable: ", a$testable, "\n", sep = "")
    cat("    Citation: ", a$citation, "\n", sep = "")
    cat("    Diagnostic: ", a$diagnostic, "\n", sep = "")
  }

  # --- Isolation read-out: deferred pending math check; see @section. ---

  invisible(x)
}


#' Format a separable_effects_assumptions Object as Markdown
#'
#' Renders the same content as `print()` but as a markdown string
#' suitable for vignette / paper reuse. Avoids drift between package
#' text and external documents.
#'
#' @param x A `"separable_effects_assumptions"` object.
#' @param style Character. Currently only `"markdown"` is supported.
#' @param ... Additional arguments (currently unused).
#' @return A single character string (one document) with markdown headers.
#' @export
format.separable_effects_assumptions <- function(x, style = "markdown", ...) {

  if (!identical(style, "markdown")) {
    stop("Only style = \"markdown\" is supported.", call. = FALSE)
  }

  lines <- character()

  push <- function(...) {
    lines <<- c(lines, paste0(..., collapse = ""))
  }

  # --- Estimand ---
  push("## Estimand")
  push("")
  push("`", x$estimand$target, "`")
  push("")
  push(x$estimand$description)
  push("")
  push("**Arms (a_Y, a_D):** ",
       paste(paste0("`", x$estimand$arms, "`"), collapse = ", "))
  push("")
  push("### ", x$estimand$decomposition_A$label)
  push("")
  push("- ", x$estimand$decomposition_A$sde)
  push("- ", x$estimand$decomposition_A$sie)
  push("")
  push("### ", x$estimand$decomposition_B$label)
  push("")
  push("- ", x$estimand$decomposition_B$sie)
  push("- ", x$estimand$decomposition_B$sde)
  push("")

  # --- Assumptions ---
  push("## Identifying assumptions")
  push("")

  emit_assumption <- function(nm, a) {
    push("### ", nm)
    push("")
    push(a$statement)
    push("")
    push("- **Testable:** ", a$testable)
    push("- **Citation:** ", a$citation)
    if (!is.null(a$diagnostic)) {
      push("- **Diagnostic:** ", a$diagnostic)
    }
    push("")
  }

  for (nm in names(x$assumptions)) {
    emit_assumption(nm, x$assumptions[[nm]])
  }
  emit_assumption("dismissible_D1", x$dismissible_D1)
  emit_assumption("dismissible_D2", x$dismissible_D2)

  # --- Isolation: deferred pending math check; see @section. ---

  paste(lines, collapse = "\n")
}
