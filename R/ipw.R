#' IPW Estimation of Separable Effects
#'
#' Estimates cumulative incidence under each treatment arm using inverse
#' probability weighting. Constructs weights that reweight observed data
#' to mimic the target intervention regime.
#'
#' @param pt_data Data frame in person-time format.
#' @param models List of fitted glm objects from [fit_hazard_models()].
#' @param treatment Character. Treatment column name.
#' @param method Character. `"ipw1"`, `"ipw2"`, or `"both"`.
#' @param cut_times Numeric vector. Time cut points.
#' @param extreme_weight_adjust Character. `"truncate"`, `"trim"`, or `"none"`.
#' @param extreme_weight_threshold Numeric. Percentile for weight adjustment.
#'
#' @return A list with:
#'   \describe{
#'     \item{cumulative_incidence}{Data frame with `k`, `arm_11`, `arm_00`,
#'       `arm_10`.}
#'     \item{weight_summary}{Data frame summarizing weight distributions.}
#'     \item{truncated_ids}{Integer vector of IDs with adjusted weights.}
#'   }
#'
#' @details
#' ## IPW Representation 1
#' Stand in A = a_Y arm. W_D reweights D-distribution across arms,
#' W_C corrects for censoring. Requires D-model + censoring model.
#'
#' ## IPW Representation 2
#' Stand in A = a_D arm. W_Y reweights Y-distribution.
#' Requires Y-model + censoring model.
#'
#' ## Weight construction
#' - `d_free_k`: interval-level D-free probability: `1 - h_D(k)`
#' - `cs_surv_d`: cumulative cause-specific D-survival: `cumprod(d_free_k)`
#'   (NOT marginal survival; conditions on being event-free)
#' - `w_d`: ratio of `cs_surv_d` across arms
#' - `w_cens`: inverse probability of remaining uncensored
#'
#' The censoring model uses observed treatment A (not overridden to a_Y),
#' because censoring weights correct for the censoring mechanism as it
#' actually operated in the observed data.
#'
#' @family internal
#' @keywords internal
ipw_estimate <- function(pt_data,
                         models,
                         treatment,
                         method,
                         cut_times,
                         extreme_weight_adjust,
                         extreme_weight_threshold) {

  id_col <- "id"
  ids <- unique(pt_data[[id_col]])
  n <- length(ids)

  # --- Predict D-free probabilities under both arms ---
  # d_free_k_ad: P(D doesn't occur at k | event-free, A_d = a_d)
  if (!is.null(models$model_d)) {
    newdata_a1 <- pt_data
    newdata_a1$A_d <- 1
    newdata_a0 <- pt_data
    newdata_a0$A_d <- 0

    pt_data$d_free_k_a1 <- 1 - stats::predict(
      models$model_d, newdata = newdata_a1, type = "response"
    )
    pt_data$d_free_k_a0 <- 1 - stats::predict(
      models$model_d, newdata = newdata_a0, type = "response"
    )

    # Cumulative cause-specific D-survival: cumprod per subject
    pt_data$cs_surv_d_a1 <- ave(
      pt_data$d_free_k_a1, pt_data[[id_col]], FUN = cumprod
    )
    pt_data$cs_surv_d_a0 <- ave(
      pt_data$d_free_k_a0, pt_data[[id_col]], FUN = cumprod
    )
  }

  # --- Predict Y-free probabilities under both arms ---
  if (!is.null(models$model_y)) {
    newdata_a1 <- pt_data
    newdata_a1$A_y <- 1
    newdata_a0 <- pt_data
    newdata_a0$A_y <- 0

    pt_data$y_free_k_a1 <- 1 - stats::predict(
      models$model_y, newdata = newdata_a1, type = "response"
    )
    pt_data$y_free_k_a0 <- 1 - stats::predict(
      models$model_y, newdata = newdata_a0, type = "response"
    )

    pt_data$cs_surv_y_a1 <- ave(
      pt_data$y_free_k_a1, pt_data[[id_col]], FUN = cumprod
    )
    pt_data$cs_surv_y_a0 <- ave(
      pt_data$y_free_k_a0, pt_data[[id_col]], FUN = cumprod
    )
  }

  # --- Censoring weights ---
  if (!is.null(models$model_c)) {
    pt_data$cens_free_k <- 1 - stats::predict(
      models$model_c, newdata = pt_data, type = "response"
    )
    pt_data$cs_surv_c <- ave(
      pt_data$cens_free_k, pt_data[[id_col]], FUN = cumprod
    )
    pt_data$w_cens <- 1 / pt_data$cs_surv_c
  }

  # --- Construct arm-specific weights ---
  # W_D: reweight D-distribution from observed arm to target arm
  if (!is.null(models$model_d)) {
    # For subjects in A=1 arm wanting D-hazard from A=0:
    # w_d = cs_surv_d_a0 / cs_surv_d_a1
    pt_data$w_d <- ifelse(
      pt_data[[treatment]] == 1,
      pt_data$cs_surv_d_a0 / pt_data$cs_surv_d_a1,
      pt_data$cs_surv_d_a1 / pt_data$cs_surv_d_a0
    )
  }

  # --- Handle extreme weights ---
  truncated_ids <- integer(0)
  if (extreme_weight_adjust != "none" && !is.null(models$model_c)) {
    w_cols <- intersect(c("w_cens", "w_d"), names(pt_data))
    for (wc in w_cols) {
      threshold_val <- stats::quantile(
        pt_data[[wc]], extreme_weight_threshold, na.rm = TRUE
      )
      extreme_rows <- which(pt_data[[wc]] > threshold_val)
      if (length(extreme_rows) > 0) {
        truncated_ids <- unique(c(
          truncated_ids,
          pt_data[[id_col]][extreme_rows]
        ))
        if (extreme_weight_adjust == "truncate") {
          pt_data[[wc]][extreme_rows] <- threshold_val
        } else if (extreme_weight_adjust == "trim") {
          pt_data <- pt_data[!pt_data[[id_col]] %in%
            pt_data[[id_col]][extreme_rows], ]
        }
      }
    }
  }

  # --- Weighted cumulative incidence ---
  # TODO: implement weighted nonparametric cumulative incidence estimator
  # using discrete_cuminc_weighted()
  cum_inc <- estimate_weighted_cum_inc(pt_data, treatment, cut_times)

  # --- Weight summary ---
  weight_summary <- summarize_weights(pt_data)

  list(
    cumulative_incidence = cum_inc,
    weight_summary = weight_summary,
    truncated_ids = truncated_ids
  )
}


#' Weighted Cumulative Incidence Estimator
#'
#' Computes cumulative incidence using weighted nonparametric estimation.
#'
#' @param pt_data Data frame with weights and event indicators.
#' @param treatment Character. Treatment column name.
#' @param cut_times Numeric vector. Time cut points.
#'
#' @return Data frame with `k`, `arm_11`, `arm_00`, `arm_10`.
#'
#' @family internal
#' @keywords internal
estimate_weighted_cum_inc <- function(pt_data, treatment, cut_times) {
  # Placeholder — will implement weighted Horvitz-Thompson / Hajek estimator
  # following the pattern from discrete_cuminc_prost() in utility_functions.R
  data.frame(
    k = cut_times,
    arm_11 = rep(NA_real_, length(cut_times)),
    arm_00 = rep(NA_real_, length(cut_times)),
    arm_10 = rep(NA_real_, length(cut_times))
  )
}


#' Summarize Weight Distributions
#'
#' @param pt_data Data frame with weight columns.
#' @return Data frame with summary statistics for each weight type.
#' @family internal
#' @keywords internal
summarize_weights <- function(pt_data) {
  w_cols <- intersect(c("w_cens", "w_d"), names(pt_data))
  if (length(w_cols) == 0) return(NULL)

  do.call(rbind, lapply(w_cols, function(wc) {
    w <- pt_data[[wc]]
    data.frame(
      weight = wc,
      mean = mean(w, na.rm = TRUE),
      median = stats::median(w, na.rm = TRUE),
      min = min(w, na.rm = TRUE),
      max = max(w, na.rm = TRUE),
      p99 = stats::quantile(w, 0.99, na.rm = TRUE),
      p999 = stats::quantile(w, 0.999, na.rm = TRUE),
      row.names = NULL
    )
  }))
}
