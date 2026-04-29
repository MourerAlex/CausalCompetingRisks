#' Classical IPW Core
#'
#' Framework-agnostic IPW math for discrete-time pooled-hazard estimators:
#' weighted hazards by interval (Hajek ratios) and the discrete-time
#' cumulative incidence formula. Vector inputs only — no dependence on
#' this package's column conventions, suitable for relocation to a future
#' shared CausalKit package.
#'
#' Used by separable IPW (Rep 1 and Rep 2) today; usable by classical
#' IPTW and any other Hajek-ratio-based hazard estimator.
#'
#' @keywords internal
NULL


#' Weighted Hazard by Interval (Hajek Ratio)
#'
#' Per-interval weighted incidence: numerator is the weighted count of
#' events at each k, denominator is the weighted count at risk at each k.
#' At-risk rows are those with non-NA `event` (the convention is that NA
#' means the subject is no longer at risk for this event type at this k).
#'
#' @param event Numeric vector. Event indicator at each row: 1 = event,
#'   0 = at risk no event, NA = no longer at risk.
#' @param k Vector of interval indices, same length as `event`.
#' @param weights Numeric vector of per-row weights, same length.
#' @return Named numeric vector: weighted hazard at each unique value of
#'   `k`, names are the k values as character.
#' @keywords internal
weighted_hazard_by_k <- function(event, k, weights) {
  at_risk <- !is.na(event)
  event   <- event[at_risk]
  k       <- k[at_risk]
  weights <- weights[at_risk]

  numer <- tapply(weights * event, k, sum, na.rm = TRUE)
  denom <- tapply(weights,         k, sum, na.rm = TRUE)
  numer / denom
}


#' Cumulative Incidence from Weighted Person-Time Data
#'
#' Discrete-time cumulative incidence of `Y` in the presence of
#' competing event `D`, computed from per-row event indicators and
#' Hajek-weighted hazards. Caller is responsible for restricting inputs
#' to a single standing-in arm.
#'
#' Computes weighted hazards by interval via [weighted_hazard_by_k()],
#' aligns to `cut_times`, and applies the standard recursion
#' \deqn{F_Y(t) = \sum_{k \le t} h_Y(k) \, (1 - h_D(k)) \, S(k-1)}
#' where \eqn{S} is the cumulative survival from joint event-freeness.
#'
#' @param y_event,d_event Numeric vectors of event indicators (Y and D).
#'   See [weighted_hazard_by_k()] for the at-risk encoding.
#' @param k Vector of interval indices.
#' @param weights Numeric vector of per-row weights.
#' @param cut_times Numeric vector of time points to evaluate at.
#' @return Numeric vector of cumulative incidence at each `cut_times[k]`.
#' @keywords internal
cum_inc_from_weighted <- function(y_event, d_event, k, weights, cut_times) {
  haz_y_by_k <- weighted_hazard_by_k(y_event, k, weights)
  haz_d_by_k <- weighted_hazard_by_k(d_event, k, weights)

  # Align to requested cut_times (missing k -> 0 hazard)
  key <- as.character(cut_times)
  haz_y_vec <- unname(haz_y_by_k[key])
  haz_d_vec <- unname(haz_d_by_k[key])
  haz_y_vec[is.na(haz_y_vec)] <- 0
  haz_d_vec[is.na(haz_d_vec)] <- 0

  # Per-interval probability of remaining event-free (joint Y and D)
  event_free <- (1 - haz_y_vec) * (1 - haz_d_vec)

  # Cumulative survival up to START of each interval (lagged)
  surv <- c(1, cumprod(event_free)[-length(event_free)])

  cumsum(haz_y_vec * (1 - haz_d_vec) * surv)
}
