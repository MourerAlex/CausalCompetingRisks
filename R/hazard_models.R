#' Fit Hazard Models
#'
#' Fits pooled logistic regression models for the Y-hazard, D-hazard, and
#' (optionally) censoring hazard on person-time data.
#'
#' @param pt_data Data frame in person-time format.
#' @param treatment Character. Treatment column name.
#' @param covariates Character vector. Covariate column names.
#' @param method Character. Estimation method (determines which models to fit).
#' @param formulas Named list or NULL. User-specified formulas.
#'
#' @return A named list of fitted glm objects:
#'   \describe{
#'     \item{model_y}{Y-hazard model (uses A_y as treatment). NULL if not needed.}
#'     \item{model_d}{D-hazard model (uses A_d as treatment). NULL if not needed.}
#'     \item{model_c}{Censoring model (uses observed A). NULL if not needed.}
#'   }
#'
#' @details
#' ## Default formula
#' `event ~ A + k + I(k^2) + I(k^3) + covariates` (additive, no interaction).
#'
#' ## Which models per method
#' - `"gformula"`: model_y + model_d
#' - `"ipw1"`: model_d + model_c
#' - `"ipw2"`: model_y + model_c
#' - `"both"`: all three
#'
#' The Y-hazard model uses `A_y` as the treatment variable; the D-hazard
#' model uses `A_d`. In observed data these are identical, but they diverge
#' in cloned datasets for cross-arm prediction.
#'
#' The censoring model always uses the observed treatment `A`, not the
#' overridden copies, because censoring weights correct for the censoring
#' mechanism as it actually operated.
#'
#' @family internal
#' @keywords internal
fit_hazard_models <- function(pt_data,
                              treatment,
                              covariates,
                              method,
                              formulas) {

  # Build default formula components
  cov_terms <- if (length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }

  time_terms <- "k + I(k^2) + I(k^3)"

  models <- list(model_y = NULL, model_d = NULL, model_c = NULL)

  # --- Y-hazard model ---
  if (method %in% c("both", "gformula", "ipw2")) {
    if (!is.null(formulas) && !is.null(formulas$y)) {
      fml_y <- formulas$y
    } else {
      rhs <- paste(c("A_y", time_terms, cov_terms), collapse = " + ")
      fml_y <- stats::as.formula(paste("y_event ~", rhs))
    }
    models$model_y <- stats::glm(
      fml_y,
      data = pt_data,
      family = stats::binomial(link = "logit")
    )
  }

  # --- D-hazard model ---
  if (method %in% c("both", "gformula", "ipw1")) {
    if (!is.null(formulas) && !is.null(formulas$d)) {
      fml_d <- formulas$d
    } else {
      rhs <- paste(c("A_d", time_terms, cov_terms), collapse = " + ")
      fml_d <- stats::as.formula(paste("d_event ~", rhs))
    }
    models$model_d <- stats::glm(
      fml_d,
      data = pt_data,
      family = stats::binomial(link = "logit")
    )
  }

  # --- Censoring model ---
  if (method %in% c("both", "ipw1", "ipw2")) {
    if (!is.null(formulas) && !is.null(formulas$c)) {
      fml_c <- formulas$c
    } else {
      # Censoring model uses observed treatment, not A_y/A_d
      rhs <- paste(c(treatment, time_terms, cov_terms), collapse = " + ")
      fml_c <- stats::as.formula(paste("c_event ~", rhs))
    }
    models$model_c <- stats::glm(
      fml_c,
      data = pt_data,
      family = stats::binomial(link = "logit")
    )
  }

  models
}
