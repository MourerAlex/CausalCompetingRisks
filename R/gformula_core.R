#' G-Formula Estimation of Separable Effects
#'
#' Estimates cumulative incidence under each treatment arm using the
#' parametric g-formula. Creates cloned datasets for arms (1,1), (0,0),
#' (1,0), and (0,1), predicts hazards using the cross-arm trick, and
#' computes per-subject cumulative incidence averaged across subjects.
#'
#' @param pt_data Data frame in person-time format. Every subject must have
#'   a row at `k = 0` (enforced by [validate_person_time()]).
#' @param models List of fitted glm objects from [fit_hazard_models()].
#' @param cut_times Numeric vector. Time cut points.
#' @param id_col Character. Subject identifier column name.
#'
#' @return A data.frame with columns `k`, `arm_11`, `arm_00`, `arm_10`,
#'   `arm_01`. The last enables Decomposition B sensitivity (see Details).
#'
#' @details
#' ## The cross-arm prediction trick
#' For arm `(a_Y, a_D)`:
#' - Predict Y-hazard using `model_y` with `A_y = a_Y`
#' - Predict D-hazard using `model_d` with `A_d = a_D`
#'
#' Arms `(1,1)` and `(0,0)` use the same treatment value in both models.
#' Arms `(1,0)` and `(0,1)` are the cross-arm regimes.
#'
#' ## Decomposition A vs B
#' - Decomposition A: direct = `arm_10 - arm_00`, indirect = `arm_11 - arm_10`
#'   (reference: `a_D = 0` in the direct)
#' - Decomposition B: direct = `arm_11 - arm_01`, indirect = `arm_01 - arm_00`
#'   (reference: `a_D = 1` in the direct)
#'
#' Both sum to the same total effect `arm_11 - arm_00`. They differ whenever
#' the `A_Y` effect on `Y` depends on `A_D` (additive `A_Y x A_D` interaction).
#'
#' ## Cumulative incidence
#' Per-subject computation, then averaged:
#' \deqn{F_Y(t) = \frac{1}{n} \sum_i \sum_{s \leq t}
#'   h_Y(s|i) (1 - h_D(s|i)) S(s-1|i)}
#'
#' where \eqn{S(s|i) = \prod_{r \leq s} (1 - h_Y(r|i))(1 - h_D(r|i))}.
#'
#' The cumprod is computed per-subject first (it does not commute with the
#' mean). The code order (mean by k, then cumsum) is mathematically equivalent
#' to (cumsum per subject, then mean) by commutativity of the double sum,
#' and is more memory-efficient.
#'
#' ## Baseline extraction
#' Uses `pt_data$k == 0` directly. Left-truncated data (subjects with no
#' `k = 0` row) is rejected upstream by [validate_person_time()] and is not
#' supported in any future version.
#'
#' @family internal
#' @keywords internal
gformula_estimate <- function(pt_data, models, cut_times, id_col) {

  # Extract baseline covariates: one row per subject at k = 0.
  # Left-truncation (subject missing k = 0) is rejected upstream.
  baseline <- pt_data[pt_data$k == 0, ]

  # --- Define arms (including arm_01 for Decomposition B sensitivity) ---
  arms <- list(
    arm_11 = c(A_y = 1, A_d = 1),
    arm_00 = c(A_y = 0, A_d = 0),
    arm_10 = c(A_y = 1, A_d = 0),
    arm_01 = c(A_y = 0, A_d = 1)
  )

  cum_inc_list <- lapply(arms, function(arm) {
    clone <- make_clone(baseline, cut_times, arm["A_y"], arm["A_d"])
    clone <- predict_hazards(clone, models)
    compute_cum_inc(clone, cut_times, id_col)
  })

  data.frame(
    k = cut_times,
    arm_11 = cum_inc_list$arm_11,
    arm_00 = cum_inc_list$arm_00,
    arm_10 = cum_inc_list$arm_10,
    arm_01 = cum_inc_list$arm_01
  )
}


#' Create a Cloned Dataset for a Treatment Arm
#'
#' Expands baseline data (one row per subject) to person-time format for a
#' specific `(a_Y, a_D)` arm configuration.
#'
#' @param baseline Data frame. One row per subject (baseline covariates).
#' @param cut_times Numeric vector. Time points.
#' @param a_y Numeric (0 or 1). Treatment value assigned to `A_y`.
#' @param a_d Numeric (0 or 1). Treatment value assigned to `A_d`.
#'
#' @return Data frame with one row per subject-time, with `A_y` and `A_d`
#'   set to the specified arm values.
#'
#' @family internal
#' @keywords internal
make_clone <- function(baseline, cut_times, a_y, a_d) {
  n <- nrow(baseline)
  n_times <- length(cut_times)

  clone <- baseline[rep(seq_len(n), each = n_times), ]
  clone$k <- rep(cut_times, times = n)
  clone$A_y <- a_y
  clone$A_d <- a_d

  clone
}


#' Predict Hazards in a Cloned Dataset
#'
#' Predicts Y-hazard and D-hazard in a cloned dataset using the fitted
#' models. The Y-model reads `A_y` and the D-model reads `A_d`, enabling
#' cross-arm prediction.
#'
#' `predict()` warnings from either model are captured and re-emitted so
#' they propagate into `fit$warnings` via `separable_effects()`'s
#' `withCallingHandlers()` wrapper. If any predicted hazards are `NA`, a
#' separate warning is emitted (silent NAs would bias cumulative incidence
#' estimates).
#'
#' @param clone Data frame. Cloned person-time data from [make_clone()].
#' @param models List of fitted glm objects.
#'
#' @return The input data frame with added columns `haz_y` and `haz_d`.
#'
#' @family internal
#' @keywords internal
predict_hazards <- function(clone, models) {
  clone$haz_y <- predict_with_warning(models$model_y, clone, "Y-hazard")
  clone$haz_d <- predict_with_warning(models$model_d, clone, "D-hazard")

  if (any(is.na(clone$haz_y))) {
    warning(
      "G-formula: Y-hazard predictions contain ", sum(is.na(clone$haz_y)),
      " NA value(s). Cumulative incidence estimates will be biased.",
      call. = FALSE
    )
  }
  if (any(is.na(clone$haz_d))) {
    warning(
      "G-formula: D-hazard predictions contain ", sum(is.na(clone$haz_d)),
      " NA value(s). Cumulative incidence estimates will be biased.",
      call. = FALSE
    )
  }

  clone
}


#' Predict with `withCallingHandlers` Warning Capture
#'
#' Internal helper wrapping `stats::predict()` so that any warnings (e.g.,
#' rank-deficient fit) are captured and re-emitted with a model label
#' prefix. This keeps the surfaces visible to `fit$warnings`.
#'
#' @param model A fitted glm object.
#' @param newdata Data frame to predict on.
#' @param label Character label prepended to re-emitted warnings.
#'
#' @return Numeric vector of predicted probabilities.
#' @family internal
#' @keywords internal
predict_with_warning <- function(model, newdata, label) {
  captured <- character()
  probs <- withCallingHandlers(
    stats::predict(model, newdata = newdata, type = "response"),
    warning = function(w) {
      captured <<- c(captured,
                     paste0(label, " predict: ", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  for (msg in captured) warning(msg, call. = FALSE)
  probs
}


#' Compute Cumulative Incidence from Predicted Hazards
#'
#' Computes the g-formula cumulative incidence of the primary event Y,
#' averaged across subjects.
#'
#' @param clone Data frame with `haz_y`, `haz_d`, `id`, `k` columns.
#' @param cut_times Numeric vector. Time points.
#'
#' @return Numeric vector of cumulative incidence at each time point.
#'
#' @details
#' At each interval k for subject i:
#' - Subdensity: `h_Y(k|i) * (1 - h_D(k|i)) * S(k-1|i)`
#' - Survival: `S(k|i) = prod_{s <= k} (1 - h_Y(s|i)) * (1 - h_D(s|i))`
#'
#' The survival (cumprod) is computed per-subject using `ave()`. The subdensity
#' is then averaged across subjects at each k, and the cumulative sum gives
#' the cumulative incidence.
#'
#' @family internal
#' @keywords internal
compute_cum_inc <- function(clone, cut_times, id_col) {

  # Per-subject event-free survival at END of each interval
  clone$surv_k <- (1 - clone$haz_y) * (1 - clone$haz_d)

  # Cumulative survival up to START of each interval (lagged cumprod)
  # S(k-1) = cumprod of surv up to previous interval; S(-1) = 1
  clone$cum_surv <- ave(
    clone$surv_k,
    clone[[id_col]],
    FUN = function(x) c(1, cumprod(x[-length(x)]))
  )

  # Subdensity at each (subject, time): h_Y * (1 - h_D) * S(k-1)
  clone$subdensity <- clone$haz_y * (1 - clone$haz_d) * clone$cum_surv

  # Average subdensity across subjects at each k, then cumulate
  mean_subdensity <- tapply(clone$subdensity, clone$k, mean)
  cumsum(mean_subdensity)
}
