#' G-Formula Estimation of Separable Effects
#'
#' Estimates cumulative incidence under each treatment arm using the
#' parametric g-formula. Creates cloned datasets for arms (1,1), (0,0),
#' and (1,0), predicts hazards using the cross-arm trick, and computes
#' per-subject cumulative incidence averaged across subjects.
#'
#' @param pt_data Data frame in person-time format.
#' @param models List of fitted glm objects from [fit_hazard_models()].
#' @param treatment Character. Treatment column name.
#' @param cut_times Numeric vector. Time cut points.
#'
#' @return A list with:
#'   \describe{
#'     \item{cumulative_incidence}{Data frame with columns `k`, `arm_11`,
#'       `arm_00`, `arm_10`.}
#'   }
#'
#' @details
#' ## The cross-arm prediction trick
#' For arm (a_Y, a_D):
#' - Predict Y-hazard using model_y with `A_y = a_Y`
#' - Predict D-hazard using model_d with `A_d = a_D`
#'
#' Arms (1,1) and (0,0) use the same treatment value for both models.
#' Arm (1,0) is the distinctive separable regime: Y-hazard under treatment,
#' D-hazard under control.
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
#' @family internal
#' @keywords internal
gformula_estimate <- function(pt_data, models, treatment, cut_times) {

  # Extract baseline data (one row per subject at k = 0)
  baseline <- pt_data[pt_data$k == min(pt_data$k), ]
  baseline <- baseline[!duplicated(baseline[["id"]]), ]

  # --- Define arms ---
  arms <- list(
    arm_11 = c(A_y = 1, A_d = 1),
    arm_00 = c(A_y = 0, A_d = 0),
    arm_10 = c(A_y = 1, A_d = 0)
  )

  cum_inc_list <- lapply(arms, function(arm) {
    clone <- make_clone(baseline, cut_times, arm["A_y"], arm["A_d"])
    clone <- predict_hazards(clone, models)
    compute_cum_inc(clone, cut_times)
  })

  data.frame(
    k = cut_times,
    arm_11 = cum_inc_list$arm_11,
    arm_00 = cum_inc_list$arm_00,
    arm_10 = cum_inc_list$arm_10
  )
}


#' Create a Cloned Dataset for a Treatment Arm
#'
#' Expands baseline data to person-time format for a specific (a_Y, a_D)
#' arm configuration.
#'
#' @param baseline Data frame. One row per subject (baseline covariates).
#' @param cut_times Numeric vector. Time points.
#' @param a_y Numeric. Treatment value for Y-hazard model (0 or 1).
#' @param a_d Numeric. Treatment value for D-hazard model (0 or 1).
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
#' @param clone Data frame. Cloned person-time data from [make_clone()].
#' @param models List of fitted glm objects.
#'
#' @return The input data frame with added columns `haz_y` and `haz_d`.
#'
#' @family internal
#' @keywords internal
predict_hazards <- function(clone, models) {
  clone$haz_y <- stats::predict(models$model_y, newdata = clone, type = "response")
  clone$haz_d <- stats::predict(models$model_d, newdata = clone, type = "response")
  clone
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
compute_cum_inc <- function(clone, cut_times) {

  # Per-subject event-free survival at END of each interval
  clone$surv_k <- (1 - clone$haz_y) * (1 - clone$haz_d)

  # Cumulative survival up to START of each interval (lagged cumprod)
  # S(k-1) = cumprod of surv up to previous interval; S(-1) = 1
  clone$cum_surv <- ave(
    clone$surv_k,
    clone[["id"]],
    FUN = function(x) c(1, cumprod(x[-length(x)]))
  )

  # Subdensity at each (subject, time): h_Y * (1 - h_D) * S(k-1)
  clone$subdensity <- clone$haz_y * (1 - clone$haz_d) * clone$cum_surv

  # Average subdensity across subjects at each k, then cumulate
  mean_subdensity <- tapply(clone$subdensity, clone$k, mean)
  cumsum(mean_subdensity)
}
