#' Inverse-Probability Weights — Core and Wrappers
#'
#' Naming convention: `ipw_*` for all weight-construction wrappers.
#' Suffix denotes the *kind* of weight (which mechanism is being inverted),
#' not the covariate set:
#'
#' * [ipw()]                  — generic core: per-row inverse of the
#'                              cumulative probability of the observed history.
#' * [ipw_cens()]             — censoring weights (Inverse Probability of
#'                              Censoring; was `ipcw`).
#' * [ipw_static_trt()]       — point-treatment weights (Inverse Probability
#'                              of Treatment, treatment fixed at baseline).
#' * [ipw_time_varying_trt()] — time-varying treatment weights — v2 placeholder.
#'
#' Stabilization is structural (two-cumprod ratio, Robins/Hernán form), not
#' a multiplicative post-process. Each wrapper accepts an optional
#' `model_num` for stabilization.
#'
#' Truncation lives downstream in [apply_weight_truncation()] and is applied
#' to the final combined IPW estimator weight (Cole & Hernán 2008 AJE), not
#' to each component.
#'
#' @keywords internal
NULL


#' Inverse-Probability Weight (Generic Core)
#'
#' Per-row inverse of the cumulative probability of the observed history:
#'
#'   W_{i,k} = prod_{j <= k} prob_num_{ij} / prod_{j <= k} prob_denom_{ij}
#'
#' With `prob_num = NULL` the numerator collapses to 1 (unstabilized form).
#' Caller is responsible for constructing the per-row probabilities; this
#' function is agnostic to which causal mechanism is being weighted.
#'
#' Cole SR, Hernán MA (2008) AJE; Hernán & Robins, *Causal Inference: What
#' If* (ch. 12, 17).
#'
#' @param prob_denom Numeric vector. P(observed history | full conditioning),
#'   one entry per row of the person-time data. Must be in (0, 1]; zeros
#'   indicate positivity violation and produce infinite weights — caller
#'   should diagnose upstream.
#' @param id Vector of subject identifiers, same length as `prob_denom`.
#' @param prob_num Optional numeric vector, same length as `prob_denom`,
#'   or NULL. P(observed history | reduced conditioning). NULL =
#'   unstabilized.
#' @param truncate Either `NULL` (default; no truncation) or a length-2
#'   numeric vector of percentile bounds. Forward-compatibility hook;
#'   v1 truncation is applied at the final combined weight downstream.
#'
#' @return Numeric weight vector, length `length(prob_denom)`.
#' @keywords internal
ipw <- function(prob_denom, id, prob_num = NULL, truncate = NULL) {
  stopifnot(
    is.numeric(prob_denom),
    length(prob_denom) == length(id),
    !any(is.na(prob_denom)),
    all(prob_denom > 0 & prob_denom <= 1)
  )
  cum_denom <- ave(prob_denom, id, FUN = cumprod)

  weight <- if (is.null(prob_num)) {
    1 / cum_denom                              # unstabilized
  } else {
    stopifnot(
      is.numeric(prob_num),
      length(prob_num) == length(prob_denom),
      !any(is.na(prob_num)),
      all(prob_num > 0 & prob_num <= 1)
    )
    cum_num <- ave(prob_num, id, FUN = cumprod)
    cum_num / cum_denom                        # stabilized
  }

  if (!is.null(truncate)) {
    qs <- stats::quantile(weight, probs = truncate, na.rm = TRUE)
    weight <- pmin(pmax(weight, qs[1]), qs[2])
  }
  weight
}


#' Inverse Probability of Censoring Weights (IPCW)
#'
#' Thin wrapper around [ipw()] for censoring weights. Predicts the censoring
#' hazard on observed covariates (no counterfactual intervention), builds
#' `prob = 1 - hazard` per row, and delegates. With `model_num` supplied,
#' returns the Robins/Hernán stabilized form (cumprod-of-numerator over
#' cumprod-of-denominator).
#'
#' Framework-agnostic: usable by classical IPTW pipelines and by separable
#' IPW pipelines for the censoring component of the weight.
#'
#' @param model_full Fitted hazard glm with full conditioning.
#' @param pt_data Person-time data frame.
#' @param id_col Character. Subject id column name.
#' @param model_num Optional fitted hazard glm with reduced conditioning
#'   (typically baseline/treatment-only). NULL = unstabilized.
#' @param truncate Forward-compat hook; default `NULL`.
#'
#' @return Numeric weight vector, length `nrow(pt_data)`.
#' @keywords internal
ipw_cens <- function(model_full, pt_data, id_col,
                     model_num = NULL, truncate = NULL) {
  haz_full   <- stats::predict(model_full, newdata = pt_data, type = "response")
  prob_denom <- 1 - haz_full
  prob_num <- if (is.null(model_num)) {
    NULL
  } else {
    1 - stats::predict(model_num, newdata = pt_data, type = "response")
  }
  ipw(prob_denom, pt_data[[id_col]], prob_num, truncate)
}


#' Inverse-Probability Weight for Static (Point) Treatment
#'
#' Builds the IPTW weight for treatment that is decided once, at baseline
#' (k = 0), and held fixed thereafter. The methodologically honest shape:
#' baseline IPTW is a per-subject scalar (one weight per individual),
#' broadcast row-wise across the subject's person-time rows.
#'
#' Predicts P(A = 1 | L) only on baseline rows (efficient: n_subjects
#' predictions, not n_rows), picks P(A_obs | L) per subject, inverts (or
#' takes the stabilized ratio with a numerator model), then broadcasts to
#' the long-format pt_data via subject id.
#'
#' Stabilization (Robins/Hernán form): provide `model_num` fit with
#' reduced conditioning, e.g. `A ~ 1` (marginal) or `A ~ V_baseline`
#' (conditional stabilization preserving effect modification on V).
#'
#' Does NOT use the cumprod-times-1 trick; baseline treatment is not
#' longitudinal and shouldn't be shoehorned through `ipw()` core.
#'
#' @param model_full Fitted glm of A on full baseline conditioning set.
#' @param pt_data Person-time data frame.
#' @param treatment_col Character. Treatment column name.
#' @param id_col Character. Subject id column name.
#' @param time_col Character. Interval index column name (default "k").
#' @param model_num Optional glm of A with reduced conditioning, for
#'   stabilization. NULL = unstabilized.
#'
#' @return Numeric weight vector, length `nrow(pt_data)`. Constant within
#'   subject (baseline treatment).
#' @keywords internal
ipw_static_trt <- function(model_full, pt_data, treatment_col, id_col,
                            time_col = "k", model_num = NULL) {
  # Predict on baseline rows only — efficient and correct for point treatment
  baseline_idx <- pt_data[[time_col]] == 0
  baseline     <- pt_data[baseline_idx, ]

  p_full     <- stats::predict(model_full, newdata = baseline, type = "response")
  p_obs_full <- ifelse(baseline[[treatment_col]] == 1, p_full, 1 - p_full)

  if (is.null(model_num)) {
    w_per_subj <- 1 / p_obs_full                          # unstabilized
  } else {
    p_n     <- stats::predict(model_num, newdata = baseline, type = "response")
    p_obs_n <- ifelse(baseline[[treatment_col]] == 1, p_n, 1 - p_n)
    w_per_subj <- p_obs_n / p_obs_full                    # stabilized ratio
  }

  # Broadcast row-wise via subject id
  names(w_per_subj) <- as.character(baseline[[id_col]])
  unname(w_per_subj[as.character(pt_data[[id_col]])])
}


#' Inverse-Probability Weight for Time-Varying Treatment (v2 — placeholder)
#'
#' Time-varying treatment IPTW: per-row P(A_k = A_obs_k | history),
#' cumulative product over time, inverted (or stabilized via a reduced-
#' conditioning numerator). Implementation deferred to v2.
#'
#' @keywords internal
ipw_time_varying_trt <- function(...) {
  stop(
    "ipw_time_varying_trt() is not implemented in v1.\n",
    "Time-varying treatment requires the Stensrud 2021 generalization. ",
    "See dev/TODO.md for the v2 plan.",
    call. = FALSE
  )
}


#' Apply Weight Truncation
#'
#' Symmetric percentile truncation (Cole & Hernán 2008 AJE) applied to
#' the source IPW weight columns on `pt_data`. Iterates over all source
#' weights actually present:
#'
#'   `c("w_cens", "w_a", "w_d_arm_10", "w_d_arm_01",`
#'    `"w_y_arm_10", "w_y_arm_01")`.
#'
#' NA values (e.g. `w_d_arm_10` on A = 0 rows) are ignored via
#' `na.rm = TRUE`. Flagged rows are recorded BEFORE clipping so the log
#' carries the raw values. A single warning summarizes how many subjects
#' and rows were affected.
#'
#' Framework-agnostic: operates on whatever weight columns are present;
#' does not know about separable arms.
#'
#' @param pt_data Person-time data frame with raw weights attached.
#' @param id_col Character. Subject id column name.
#' @param truncate Either NULL (no truncation) or a length-2 numeric
#'   vector of percentile bounds, e.g. `c(0.01, 0.99)`. Setting the
#'   lower bound to 0 reduces to upper-tail-only truncation.
#'
#' @return A list with:
#'   \describe{
#'     \item{pt_data}{Adjusted person-time data with weights clipped.}
#'     \item{flagged_ids}{Vector of unique subject IDs whose weights
#'       exceeded the bounds across any weight column.}
#'     \item{flagged_log}{Long-format data.frame with one row per
#'       flagged `(subject, interval, weight column)` triple.}
#'   }
#' @keywords internal
apply_weight_truncation <- function(pt_data, id_col, truncate = NULL) {

  empty_log <- data.frame(
    id     = pt_data[[id_col]][0],
    weight = character(0),
    k      = numeric(0),
    value  = numeric(0),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (is.null(truncate)) {
    return(list(
      pt_data     = pt_data,
      flagged_ids = integer(0),
      flagged_log = empty_log
    ))
  }

  w_cols <- intersect(
    c("w_cens", "w_a",
      "w_d_arm_10", "w_d_arm_01",
      "w_y_arm_10", "w_y_arm_01"),
    names(pt_data)
  )

  log_rows <- list()

  for (wc in w_cols) {
    qs <- stats::quantile(pt_data[[wc]], probs = truncate, na.rm = TRUE)
    lower_val <- qs[1]
    upper_val <- qs[2]

    extreme_rows <- which(
      pt_data[[wc]] < lower_val | pt_data[[wc]] > upper_val
    )

    if (length(extreme_rows) > 0) {
      log_rows[[wc]] <- data.frame(
        id     = pt_data[[id_col]][extreme_rows],
        weight = wc,
        k      = pt_data$k[extreme_rows],
        value  = pt_data[[wc]][extreme_rows],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      # Symmetric clip (Cole & Hernán 2008)
      pt_data[[wc]] <- pmin(pmax(pt_data[[wc]], lower_val), upper_val)
    }
  }

  flagged_log <- if (length(log_rows) == 0) empty_log
                 else do.call(rbind, log_rows)
  flagged_ids   <- unique(flagged_log$id)
  n_flagged_rows <- nrow(flagged_log)

  if (length(flagged_ids) > 0) {
    warning(
      length(flagged_ids), " subject(s) had weights truncated at ",
      "p=[", truncate[1], ", ", truncate[2], "] (",
      n_flagged_rows, " row(s) affected). ",
      "See $flagged_log for details.",
      call. = FALSE
    )
  }

  list(
    pt_data     = pt_data,
    flagged_ids = flagged_ids,
    flagged_log = flagged_log
  )
}


#' Summarize Weight Distributions
#'
#' Reports distributional statistics for every source weight column on
#' `pt_data`, for both raw (`*_raw`) and truncated versions when both
#' exist. Source-weight base names recognized:
#'
#'   `c("w_cens", "w_a", "w_d_arm_10", "w_d_arm_01",`
#'    `"w_y_arm_10", "w_y_arm_01")`.
#'
#' @param pt_data Data frame with weight columns.
#' @return Data frame with one row per weight column (raw and truncated
#'   for each base name when present) and columns
#'   `weight, n_nonNA, mean, median, min, p99, p999, max`. NULL if no
#'   weight columns are present.
#' @keywords internal
summarize_weights <- function(pt_data) {
  base_names <- c("w_cens", "w_a",
                  "w_d_arm_10", "w_d_arm_01",
                  "w_y_arm_10", "w_y_arm_01")

  w_cols <- character()
  for (bn in base_names) {
    if (paste0(bn, "_raw") %in% names(pt_data)) {
      w_cols <- c(w_cols, paste0(bn, "_raw"))
    }
    if (bn %in% names(pt_data)) {
      w_cols <- c(w_cols, bn)
    }
  }

  if (length(w_cols) == 0) return(NULL)

  do.call(rbind, lapply(w_cols, function(wc) {
    w <- pt_data[[wc]]
    data.frame(
      weight  = wc,
      n_nonNA = sum(!is.na(w)),
      mean    = mean(w, na.rm = TRUE),
      median  = stats::median(w, na.rm = TRUE),
      min     = min(w, na.rm = TRUE),
      p99     = stats::quantile(w, 0.99,  na.rm = TRUE),
      p999    = stats::quantile(w, 0.999, na.rm = TRUE),
      max     = max(w, na.rm = TRUE),
      row.names = NULL
    )
  }))
}
