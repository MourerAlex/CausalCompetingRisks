#' Fit Hazard Models
#'
#' Fits pooled logistic regression models for the Y-hazard, D-hazard, and
#' (optionally) censoring hazard on person-time data. For each fitted model,
#' collects diagnostics (convergence, min/max fitted probability, positivity
#' flag, captured `glm()` warnings).
#'
#' @param pt_data Data frame in person-time format.
#' @param treatment Character. Treatment column name (used for the censoring
#'   model; Y and D models use `A_y` / `A_d` working copies).
#' @param covariates Character vector. Covariate column names.
#' @param active_methods Character vector. Subset of
#'   `c("gformula", "ipw1", "ipw2")` indicating which methods will run.
#'   Determines which models get fit.
#' @param formulas Named list or NULL. User-specified formulas (names `y`,
#'   `d`, `c`). Any entry absent falls back to the default formula.
#' @param censoring_weights Logical. When FALSE, the censoring model is not
#'   fit and `model_c` stays NULL.
#'
#' @return A named list with two entries:
#'   \describe{
#'     \item{models}{Named list: `model_y`, `model_d`, `model_c` (glm
#'       objects or NULL).}
#'     \item{checks}{Named list: `y`, `d`, `c` (per-model diagnostics or
#'       NULL). See [check_fitted_positivity()] for the per-model structure.}
#'   }
#'
#' @details
#' ## Default formula
#' `flag ~ A_y/A_d/treatment + k + I(k^2) + I(k^3) + covariates` (additive,
#' no interaction).
#'
#' ## Which models per method
#' - `"gformula"`: model_y + model_d
#' - `"ipw1"`: model_d + model_c
#' - `"ipw2"`: model_y + model_c
#'
#' ## Treatment handling
#' - Y-hazard model: uses `A_y` (a working-copy column on pt_data)
#' - D-hazard model: uses `A_d` (same idea)
#' - Censoring model: uses the observed `treatment` column directly
#'
#' In observed data all three are identical; `A_y`/`A_d` diverge only in
#' cloned datasets used for cross-arm prediction downstream.
#'
#' @family internal
#' @keywords internal
fit_hazard_models <- function(pt_data,
                              treatment,
                              covariates,
                              active_methods,
                              formulas,
                              censoring_weights = TRUE) {

  # Build default formula components
  cov_terms <- if (length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }

  time_terms <- "k + I(k^2) + I(k^3)"

  models <- list(model_y = NULL, model_d = NULL, model_c = NULL)
  checks <- list(y = NULL, d = NULL, c = NULL)

  # --- Y-hazard model (needed for gformula and ipw2) ---
  if (any(c("gformula", "ipw2") %in% active_methods)) {
    fml_y <- formulas$y %||% stats::as.formula(
      paste("y_flag ~", paste(c("A_y", time_terms, cov_terms),
                              collapse = " + "))
    )
    fit_result <- fit_logistic(fml_y, pt_data, "Y-hazard")
    models$model_y <- fit_result$model
    checks$y <- fit_result$check
  }

  # --- D-hazard model (needed for gformula and ipw1) ---
  if (any(c("gformula", "ipw1") %in% active_methods)) {
    fml_d <- formulas$d %||% stats::as.formula(
      paste("d_flag ~", paste(c("A_d", time_terms, cov_terms),
                              collapse = " + "))
    )
    fit_result <- fit_logistic(fml_d, pt_data, "D-hazard")
    models$model_d <- fit_result$model
    checks$d <- fit_result$check
  }

  # --- Censoring model (only if IPW requested AND censoring_weights = TRUE) ---
  if (any(c("ipw1", "ipw2") %in% active_methods) && censoring_weights) {
    fml_c <- formulas$c %||% stats::as.formula(
      paste("c_flag ~", paste(c(treatment, time_terms, cov_terms),
                              collapse = " + "))
    )
    fit_result <- fit_logistic(fml_c, pt_data, "C-hazard")
    models$model_c <- fit_result$model
    checks$c <- fit_result$check
  }

  list(models = models, checks = checks)
}


#' Fit a Logistic GLM with Warning Capture and Positivity Check
#'
#' Internal helper. Fits a binomial logistic glm, captures any `glm()`
#' warnings via `withCallingHandlers()`, re-emits them so `separable_effects()`'s
#' collector grabs them, and builds a diagnostics list via
#' [check_fitted_positivity()].
#'
#' Used by [fit_hazard_models()] (Y, D, C hazards on person-time) and by
#' [fit_propensity()] (treatment model on baseline rows). Family-agnostic
#' callers can extend the `family` arg later if needed.
#'
#' @param formula A fitted model formula.
#' @param data Data frame passed to `glm()`.
#' @param label Human-readable label for warnings and diagnostics.
#'
#' @return A list with `model` (the glm object) and `check` (diagnostics).
#' @family internal
#' @keywords internal
fit_logistic <- function(formula, data, label) {
  glm_warnings <- character()

  model <- withCallingHandlers(
    stats::glm(formula, data = data,
               family = stats::binomial(link = "logit")),
    warning = function(w) {
      glm_warnings <<- c(glm_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # Re-emit glm warnings so they are collected by separable_effects()'s handler
  for (w in glm_warnings) warning(w, call. = FALSE)

  check <- check_fitted_positivity(model, label, glm_warnings = glm_warnings)

  list(model = model, check = check)
}


#' Check Fitted Probabilities (Positivity Signal)
#'
#' Inspects predicted probabilities from a fitted glm and returns a
#' diagnostics list with the continuous signal (min and max fitted). Emits
#' warnings for: non-convergence, NA fitted values, and predicted
#' probabilities near 0 or 1 (the user interprets the continuous extremes;
#' no binary positivity-violation flag is reported because any specific
#' threshold would be arbitrary).
#'
#' @param model A fitted glm object.
#' @param model_label Character label for the model (e.g., "C-hazard").
#' @param warn_eps Threshold below 0 (and above 1-warn_eps) for emitting an
#'   "extreme probabilities predicted" warning. Default 1e-6. Only controls
#'   whether a warning fires — the continuous min/max are always returned.
#' @param glm_warnings Character vector of glm warnings captured at fit time.
#'
#' @return A list with:
#'   \describe{
#'     \item{label}{Model label.}
#'     \item{converged}{Logical.}
#'     \item{min_fitted, max_fitted}{Extremes of fitted probabilities (or NA
#'       if fitted values themselves contain NAs). Always reported so users
#'       can judge positivity for themselves.}
#'     \item{glm_warnings}{Character vector passed through from the fit.}
#'   }
#' @family internal
#' @keywords internal
check_fitted_positivity <- function(model, model_label,
                                    warn_eps = 1e-6,
                                    glm_warnings = character()) {

  converged <- isTRUE(model$converged)
  if (!converged) {
    warning(
      model_label, " model did not converge.",
      call. = FALSE
    )
  }

  probs <- stats::fitted(model)

  if (any(is.na(probs))) {
    warning(
      model_label, " model has NA fitted values - fit may have failed.",
      call. = FALSE
    )
    return(list(
      label = model_label,
      converged = converged,
      min_fitted = NA_real_,
      max_fitted = NA_real_,
      glm_warnings = glm_warnings
    ))
  }

  min_p <- min(probs)
  max_p <- max(probs)

  # Warn about extreme probabilities (user interprets what "extreme" means).
  # No binary positivity-violation flag is stored — the min/max are the
  # continuous signal, and any cutoff would be arbitrary.
  if (min_p < warn_eps || max_p > 1 - warn_eps) {
    warning(
      model_label,
      " model predicted extreme probabilities ",
      "(min=", signif(min_p, 3), ", max=", signif(max_p, 3), "). ",
      "Inspect diagnostic(fit) to assess positivity.",
      call. = FALSE
    )
  }

  list(
    label = model_label,
    converged = converged,
    min_fitted = min_p,
    max_fitted = max_p,
    glm_warnings = glm_warnings
  )
}


# ==============================================================================
# Hazard prediction & survival primitives (framework-agnostic)
# ==============================================================================

#' Predict Hazard Under a Counterfactual Treatment Assignment
#'
#' Replaces `treatment_var` in `data` with `treatment_value`, then calls
#' `stats::predict()` on `model`. Returns predicted hazards aligned to rows.
#' Returns NULL if `model` is NULL. Knows nothing about Y, D, separable arms,
#' or any specific causal framework.
#'
#' @param model Fitted hazard glm (binomial), or NULL.
#' @param data Person-time data frame.
#' @param treatment_var Character. Name of the column to overwrite.
#' @param treatment_value Numeric. Counterfactual value to assign.
#' @return Numeric vector of predicted hazards, or NULL.
#' @keywords internal
predict_hazard_under <- function(model, data, treatment_var, treatment_value) {
  if (is.null(model)) return(NULL)
  newdata <- data
  newdata[[treatment_var]] <- treatment_value
  stats::predict(model, newdata = newdata, type = "response")
}


#' Per-Subject Cumulative Product of Survival Probability
#'
#' Computes prod_{j <= s}(1 - haz_j) within each subject id. Returned vector
#' is aligned to the input rows. A core survival-analysis primitive.
#'
#' @param haz Numeric vector of hazards.
#' @param id Subject id vector (same length as `haz`).
#' @return Numeric vector of cumulative survival, aligned to input rows.
#' @keywords internal
cumprod_survival <- function(haz, id) {
  ave(1 - haz, id, FUN = cumprod)
}
