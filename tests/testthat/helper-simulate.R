# ============================================================================
# Test helper: simulate_separable_data()
#
# Discrete-time DGP for the separable-effects test ladder. Generates
# subject-level competing-events data with known per-arm cumulative
# incidence (analytical truth) so tests can verify estimator behaviour
# end-to-end.
#
# Covariate model (option A — locked 2026-04-29 with user):
#   L_j ~ Bernoulli(prev_L_j),  independent across j = 1..d
#
# Treatment model:
#   A | L ~ Bernoulli(plogis(beta_a$intercept + beta_a$L · L))
#
# In observed data, the two component-treatments are glued: A_Y = A_D = A.
# The simulator records `A_Y` and `A_D` as separate columns to make the
# decomposition explicit; tests that toggle Δ1 / Δ2 control the dependence
# of the hazard models on these via beta_y$A_D and beta_d$A_Y.
#
# Hazard models (discrete time, k = 0, 1, ..., K-1):
#   h_Y(k | L, A_Y, A_D) = plogis(b$intercept + b$L · L + b$A_Y · A_Y + b$A_D · A_D + b$k · k)
#   h_D(k | L, A_Y, A_D) = same shape with beta_d
#   h_C(k | L, A)        = plogis(c$intercept + c$L · L + c$A · A + c$k · k)
#
# Temporal ordering within an interval: C first, then D, then Y.
# Subjects who reach k = K-1 without an event are administratively
# censored (event_type = c, event_time = K).
#
# Returns:
#   list(
#     data   = data.frame(id, event_time, event_type, A, L_1, ..., L_d),
#     truth  = list(arm_cif, contrasts),    # arm CIFs + 5 contrasts at each k
#     params = list(...)
#   )
# ============================================================================

#' Simulate Subject-Level Competing-Events Data with Known Truth
#'
#' @param n Integer. Number of subjects.
#' @param K Integer. Number of intervals (k = 0, 1, ..., K-1; admin censoring
#'   at K).
#' @param prev_L Numeric vector of length d. Marginal P(L_j = 1).
#' @param beta_y Named list of Y-hazard coefficients: intercept (numeric),
#'   L (length d), A_Y, A_D, k.
#' @param beta_d Named list of D-hazard coefficients: same shape as beta_y.
#' @param beta_c Named list of C-hazard coefficients: intercept, L (length d),
#'   A, k.
#' @param beta_a Named list of treatment-propensity coefficients: intercept,
#'   L (length d).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return list(data, truth, params). See file header.
simulate_separable_data <- function(n,
                                    K       = 12,
                                    prev_L  = c(0.5),
                                    beta_y  = list(intercept = -3,   L = c(0),  A_Y = 0, A_D = 0, k = 0),
                                    beta_d  = list(intercept = -3.5, L = c(0),  A_Y = 0, A_D = 0, k = 0),
                                    beta_c  = list(intercept = -4,   L = c(0),  A   = 0,         k = 0),
                                    beta_a  = list(intercept = 0,    L = c(0)),
                                    seed    = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # --- Input validation ---
  d <- length(prev_L)
  stopifnot(d >= 0, n >= 1, K >= 1)
  stopifnot(length(beta_y$L) == d, length(beta_d$L) == d,
            length(beta_c$L) == d, length(beta_a$L) == d)
  for (b_name in c("beta_y", "beta_d")) {
    b <- get(b_name)
    stopifnot(all(c("intercept", "L", "A_Y", "A_D", "k") %in% names(b)))
  }
  stopifnot(all(c("intercept", "L", "A", "k") %in% names(beta_c)))
  stopifnot(all(c("intercept", "L") %in% names(beta_a)))

  # --- Draw covariates and treatment ---
  L_mat <- matrix(NA_integer_, nrow = n, ncol = d)
  if (d > 0) {
    for (j in seq_len(d)) {
      L_mat[, j] <- stats::rbinom(n, 1, prev_L[j])
    }
    colnames(L_mat) <- paste0("L_", seq_len(d))
  }

  if (d > 0) {
    pi_a <- stats::plogis(beta_a$intercept + as.numeric(L_mat %*% beta_a$L))
  } else {
    pi_a <- rep(stats::plogis(beta_a$intercept), n)
  }
  A <- stats::rbinom(n, 1, pi_a)

  # --- Per-subject discrete-time event simulation ---
  event_time <- integer(n)
  event_type <- integer(n)  # 0 = c (censored), 1 = y, 2 = d

  hazard_logit <- function(b, L_row, a_y, a_d, k) {
    eta <- b$intercept + b$A_Y * a_y + b$A_D * a_d + b$k * k
    if (d > 0) eta <- eta + sum(b$L * L_row)
    stats::plogis(eta)
  }
  hazard_logit_c <- function(L_row, a, k) {
    eta <- beta_c$intercept + beta_c$A * a + beta_c$k * k
    if (d > 0) eta <- eta + sum(beta_c$L * L_row)
    stats::plogis(eta)
  }

  for (i in seq_len(n)) {
    L_row  <- if (d > 0) L_mat[i, ] else integer(0)
    a_obs  <- A[i]
    # In observed data, A_Y = A_D = A:
    a_y    <- a_obs
    a_d    <- a_obs

    finished <- FALSE
    for (k in 0:(K - 1L)) {
      h_c <- hazard_logit_c(L_row, a_obs, k)
      h_d <- hazard_logit(beta_d, L_row, a_y, a_d, k)
      h_y <- hazard_logit(beta_y, L_row, a_y, a_d, k)

      # Temporal ordering: C, then D, then Y
      if (stats::runif(1) < h_c) {
        event_time[i] <- k; event_type[i] <- 0L; finished <- TRUE; break
      }
      if (stats::runif(1) < h_d) {
        event_time[i] <- k; event_type[i] <- 2L; finished <- TRUE; break
      }
      if (stats::runif(1) < h_y) {
        event_time[i] <- k; event_type[i] <- 1L; finished <- TRUE; break
      }
    }
    if (!finished) {
      event_time[i] <- K
      event_type[i] <- 0L  # admin censoring at K
    }
  }

  data_df <- data.frame(
    id         = seq_len(n),
    event_time = event_time,
    event_type = event_type,
    A          = A,
    stringsAsFactors = FALSE
  )
  if (d > 0) {
    data_df <- cbind(data_df, as.data.frame(L_mat))
  }

  # --- Analytical truth: marginal arm CIFs ---
  truth <- compute_arm_cif_truth(
    K = K, prev_L = prev_L,
    beta_y = beta_y, beta_d = beta_d
  )

  list(
    data   = data_df,
    truth  = truth,
    params = list(n = n, K = K, prev_L = prev_L,
                  beta_y = beta_y, beta_d = beta_d,
                  beta_c = beta_c, beta_a = beta_a)
  )
}


#' Closed-Form Marginal Arm CIFs and Contrasts
#'
#' Enumerates 2^d binary L-strata, computes conditional CIF per stratum
#' and arm via the discrete recursion, weights by stratum probability,
#' and assembles the marginal arm CIFs and the five contrasts (total,
#' SDE-A, SIE-A, SDE-B, SIE-B) at each k.
#'
#' @param K Number of intervals.
#' @param prev_L Marginal Bernoulli probabilities of L (length d).
#' @param beta_y,beta_d Hazard coefficient lists (see simulate_separable_data).
#' @return list(arm_cif, contrasts).
compute_arm_cif_truth <- function(K, prev_L, beta_y, beta_d) {

  d     <- length(prev_L)
  arms  <- list(arm_11 = c(1, 1), arm_00 = c(0, 0),
                arm_10 = c(1, 0), arm_01 = c(0, 1))
  k_seq <- 0:(K - 1L)

  # Enumerate strata as a (2^d x d) integer matrix
  if (d == 0) {
    L_strata <- matrix(integer(0), nrow = 1, ncol = 0)
    p_stratum <- 1
  } else {
    L_strata <- as.matrix(expand.grid(rep(list(c(0L, 1L)), d)))
    colnames(L_strata) <- paste0("L_", seq_len(d))
    p_stratum <- apply(L_strata, 1, function(row) {
      prod(prev_L^row * (1 - prev_L)^(1 - row))
    })
  }

  haz <- function(b, L_row, a_y, a_d, k) {
    eta <- b$intercept + b$A_Y * a_y + b$A_D * a_d + b$k * k
    if (d > 0) eta <- eta + sum(b$L * L_row)
    stats::plogis(eta)
  }

  # CIF at each k for one stratum × one arm: F_Y(k) = cumsum_{j<=k} h_Y(j) (1-h_D(j)) S(j-1)
  cif_one <- function(L_row, a_y, a_d) {
    h_y <- vapply(k_seq, function(k) haz(beta_y, L_row, a_y, a_d, k),
                  numeric(1))
    h_d <- vapply(k_seq, function(k) haz(beta_d, L_row, a_y, a_d, k),
                  numeric(1))
    # S(k) = prod_{j<=k} (1 - h_y(j)) (1 - h_d(j)); S(-1) = 1
    one_minus_step <- (1 - h_y) * (1 - h_d)
    # lag-cumprod: at index k, value = prod up to (k-1)
    lag_cumprod <- c(1, cumprod(one_minus_step)[-K])
    # Cause-specific subdensity for Y at k: h_y(k) * (1 - h_d(k)) * S(k-1)
    contrib <- h_y * (1 - h_d) * lag_cumprod
    cumsum(contrib)
  }

  arm_cif_list <- lapply(arms, function(arm_vec) {
    cif_strata <- vapply(seq_len(nrow(L_strata)), function(s) {
      L_row <- if (d > 0) L_strata[s, ] else integer(0)
      cif_one(L_row, a_y = arm_vec[1], a_d = arm_vec[2])
    }, numeric(K))
    # cif_strata is K x n_strata — marginalise
    if (d == 0) {
      as.numeric(cif_strata)
    } else {
      as.numeric(cif_strata %*% p_stratum)
    }
  })

  arm_cif <- data.frame(
    k      = k_seq,
    arm_11 = arm_cif_list$arm_11,
    arm_00 = arm_cif_list$arm_00,
    arm_10 = arm_cif_list$arm_10,
    arm_01 = arm_cif_list$arm_01
  )

  contrasts <- data.frame(
    k     = k_seq,
    total = arm_cif$arm_11 - arm_cif$arm_00,
    sde_a = arm_cif$arm_10 - arm_cif$arm_00,
    sie_a = arm_cif$arm_11 - arm_cif$arm_10,
    sde_b = arm_cif$arm_11 - arm_cif$arm_01,
    sie_b = arm_cif$arm_01 - arm_cif$arm_00
  )

  list(arm_cif = arm_cif, contrasts = contrasts)
}
