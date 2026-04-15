#' Estimate Separable Effects for Competing Events
#'
#' Main entry point for the package. Fits discrete-time pooled logistic
#' regression models and estimates total, separable direct, and separable
#' indirect effects on the cumulative incidence of a primary event in the
#' presence of a competing event.
#'
#' @param data A data.frame. Can be subject-level (one row per subject) or
#'   person-time (one row per subject-interval). Detected automatically.
#' @param id Character. Name of the subject identifier column.
#' @param time Character. Name of the event/censoring time column.
#' @param event Character. Name of the event type column. Must have exactly
#'   3 unique values (censoring, primary event Y, competing event D).
#' @param treatment Character. Name of the binary treatment column (0/1).
#' @param covariates Character vector. Names of baseline covariate columns.
#' @param event_y Character or NULL. Value in `event` indicating the primary
#'   event. If NULL, uses integer convention: 0 = censored, 1 = Y, 2 = D.
#' @param event_d Character or NULL. Value indicating the competing event.
#' @param event_c Character or NULL. Value indicating censoring.
#' @param method Character. Estimation method: `"both"` (default), `"gformula"`,
#'   `"ipw1"`, or `"ipw2"`.
#' @param times Numeric. Time grid for person-time expansion. NULL = 12
#'   equally-spaced intervals; scalar = that many intervals; vector = cut points.
#' @param eval_times Numeric vector. Time points at which to report contrasts.
#'   NULL = all time points.
#' @param formulas Named list or NULL. Override default model formulas. Names
#'   should be `"y"`, `"d"`, and/or `"c"`.
#' @param time_varying Character vector or NULL. Not implemented in v1; raises
#'   an error if non-NULL.
#' @param extreme_weight_adjust Character. Weight adjustment method for IPW:
#'   `"truncate"` (default), `"trim"`, or `"none"`.
#' @param extreme_weight_threshold Numeric. Percentile threshold for weight
#'   adjustment (default 0.999).
#' @param bootstrap Logical. Whether to compute bootstrap confidence intervals.
#' @param n_boot Integer. Number of bootstrap replicates (default 500).
#' @param alpha Numeric. Significance level for confidence intervals
#'   (default 0.05).
#'
#' @return An S3 object of class `"causal_cr"` containing:
#'   \describe{
#'     \item{cumulative_incidence}{Data frame with cumulative incidence for
#'       each arm (1,1), (0,0), (1,0) at each time point.}
#'     \item{contrasts}{Data frame with risk differences and risk ratios for
#'       total, separable direct, and separable indirect effects.}
#'     \item{models}{List of fitted glm objects.}
#'     \item{person_time}{The expanded person-time dataset.}
#'     \item{diagnostics}{Weight summaries and isolation checks.}
#'     \item{bootstrap}{NULL or list with bootstrap replicates and CIs.}
#'   }
#'
#' @references
#' Stensrud MJ, Young JG, Didelez V, Robins JM, Hernan MA (2020).
#' "Separable Effects for Causal Inference in the Presence of Competing Events."
#' *Journal of the American Statistical Association*, 115(531), 1363-1375.
#' \doi{10.1080/01621459.2020.1765783}
#'
#' Stensrud MJ, Hernan MA, Tchetgen Tchetgen EJ, Robins JM, Didelez V,
#' Young JG (2021). "A Generalized Theory of Separable Effects in Competing
#' Event Settings." *Lifetime Data Analysis*, 27(4), 588-631.
#' \doi{10.1007/s10985-021-09530-8}
#'
#' @family main
#' @seealso [risk()], [contrast()], [diagnostic()]
#'
#' @examples
#' \dontrun{
#' data(prostate_data)
#' fit <- causal_cr(
#'   data = prostate_data,
#'   id = "id",
#'   time = "event_time",
#'   event = "event_type",
#'   treatment = "A",
#'   covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin")
#' )
#' summary(fit)
#' plot(risk(fit))
#' plot(contrast(fit))
#' }
#'
#' @export
causal_cr <- function(data,
                      id = "id",
                      time = "time",
                      event = "event_type",
                      treatment = "A",
                      covariates = character(),
                      event_y = NULL,
                      event_d = NULL,
                      event_c = NULL,
                      method = "both",
                      times = NULL,
                      eval_times = NULL,
                      formulas = NULL,
                      time_varying = NULL,
                      extreme_weight_adjust = "truncate",
                      extreme_weight_threshold = 0.999,
                      bootstrap = FALSE,
                      n_boot = 500L,
                      alpha = 0.05) {

  # Capture call
  cl <- match.call()

  # --- Input validation ---
  validate_input(
    data = data,
    id = id,
    time = time,
    event = event,
    treatment = treatment,
    covariates = covariates,
    event_y = event_y,
    event_d = event_d,
    event_c = event_c,
    method = method,
    time_varying = time_varying
  )

  # --- Data preparation ---
  prep <- to_person_time(
    data = data,
    id = id,
    time = time,
    event = event,
    treatment = treatment,
    covariates = covariates,
    event_y = event_y,
    event_d = event_d,
    event_c = event_c,
    times = times
  )

  pt_data <- prep$person_time
  cut_times <- prep$cut_times

  # --- Fit hazard models ---
  models <- fit_hazard_models(
    pt_data = pt_data,
    treatment = treatment,
    covariates = covariates,
    method = method,
    formulas = formulas
  )

  # --- Estimation ---
  results <- list()

  if (method %in% c("both", "gformula")) {
    results$gformula <- gformula_estimate(
      pt_data = pt_data,
      models = models,
      treatment = treatment,
      cut_times = cut_times
    )
  }

  if (method %in% c("both", "ipw1", "ipw2")) {
    results$ipw <- ipw_estimate(
      pt_data = pt_data,
      models = models,
      treatment = treatment,
      method = method,
      cut_times = cut_times,
      extreme_weight_adjust = extreme_weight_adjust,
      extreme_weight_threshold = extreme_weight_threshold
    )
  }

  # --- Contrasts ---
  cumulative_incidence <- build_cumulative_incidence(results, method)
  contrasts <- compute_contrasts(cumulative_incidence, eval_times)

  # --- Diagnostics ---
  diagnostics <- build_diagnostics(results, method)

  # --- Bootstrap ---
  boot_results <- NULL
  if (bootstrap) {
    boot_results <- bootstrap_causal_cr(
      data = data,
      id = id,
      time = time,
      event = event,
      treatment = treatment,
      covariates = covariates,
      event_y = event_y,
      event_d = event_d,
      event_c = event_c,
      method = method,
      times = times,
      eval_times = eval_times,
      formulas = formulas,
      extreme_weight_adjust = extreme_weight_adjust,
      extreme_weight_threshold = extreme_weight_threshold,
      n_boot = n_boot,
      alpha = alpha
    )
  }

  # --- Build return object ---
  structure(
    list(
      cumulative_incidence = cumulative_incidence,
      contrasts = contrasts,
      models = models,
      person_time = pt_data,
      diagnostics = diagnostics,
      bootstrap = boot_results,
      call = cl,
      method = method,
      times = cut_times,
      n = length(unique(pt_data[[id]]))
    ),
    class = "causal_cr"
  )
}
