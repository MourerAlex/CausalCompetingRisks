# ============================================================================
# Steps 1-3 of the user-locked DGP test ladder
# (see memory: project_test_ladder.md).
#
# Step 1 — pooled logistic recovers true Y / D / C hazard coefficients on
#          a confounded DGP with one binary covariate
# Step 2 — discrete S(k) = prod(1 - h_y(j)) matches KM with NO censoring
# Step 3 — same with INDEPENDENT censoring
#
# Steps 1-3 test primitives below the package's main entry. They fit
# glms directly on hand-built person-time data so we can compare against
# survival::survfit at the level of cumulative survival, which is what
# the package's CIF ultimately reduces to in degenerate cases.
# ============================================================================

# ---- Local helper: build person-time manually ------------------------------
# Bypasses to_person_time() (which requires exactly 3 event types in the
# data) so that Y-only and Y+C-only DGPs can be tested cleanly.
build_pt <- function(data_df, K) {
  rows <- vector("list", nrow(data_df))
  for (i in seq_len(nrow(data_df))) {
    et  <- data_df$event_time[i]
    ety <- data_df$event_type[i]
    final_k <- min(et, K - 1L)
    k_seq <- 0:final_k
    is_terminal <- (k_seq == et) & (et <= (K - 1L))
    rows[[i]] <- data.frame(
      id     = data_df$id[i],
      k      = k_seq,
      y_flag = as.integer(is_terminal & ety == 1L),
      d_flag = as.integer(is_terminal & ety == 2L),
      c_flag = as.integer(is_terminal & ety == 0L & et < K),
      stringsAsFactors = FALSE
    )
    # Apply the temporal-ordering NA convention used in the package
    rows[[i]]$d_flag[rows[[i]]$c_flag == 1L] <- NA_integer_
    rows[[i]]$y_flag[is.na(rows[[i]]$d_flag) | rows[[i]]$d_flag == 1L] <- NA_integer_
  }
  out <- do.call(rbind, rows)
  baseline_cols <- setdiff(names(data_df), c("event_time", "event_type"))
  out <- merge(out, data_df[, baseline_cols, drop = FALSE], by = "id",
               sort = FALSE)
  out[order(out$id, out$k), ]
}


# ============================================================================
# Step 1: pooled logistic recovers true hazard coefficients
# ============================================================================
# In observed data A_Y = A_D = A. So glm fits one A coefficient per hazard.
# Under Δ1 (beta_y$A_D = 0) the Y-coefficient on A recovers beta_y$A_Y.
# Under Δ2 (beta_d$A_Y = 0) the D-coefficient on A recovers beta_d$A_D.
# We set both isolation conditions in this DGP so each true coefficient is
# uniquely identified. n is large to keep MLE noise small.
# ============================================================================

test_that("Step 1: pooled logistic recovers true Y / D / C hazard coefficients", {
  skip_on_cran()

  beta_y_true <- list(intercept = -3,   L = c(0.4),  A_Y = -0.6, A_D = 0,    k = -0.05)
  beta_d_true <- list(intercept = -3.5, L = c(0.3),  A_Y = 0,    A_D = -0.4, k = 0)
  beta_c_true <- list(intercept = -4,   L = c(0.2),  A   = 0.1,              k = 0)

  sim <- simulate_separable_data(
    n      = 8000, K = 12, prev_L = c(0.5),
    beta_y = beta_y_true,
    beta_d = beta_d_true,
    beta_c = beta_c_true,
    beta_a = list(intercept = 0, L = c(0)),     # no confounding for step 1
    seed   = 1
  )

  pt <- build_pt(sim$data, K = 12)
  pt$A_y <- pt$A
  pt$A_d <- pt$A

  fit_y <- stats::glm(y_flag ~ A_y + L_1 + k,
                      data = pt, family = stats::binomial())
  fit_d <- stats::glm(d_flag ~ A_d + L_1 + k,
                      data = pt, family = stats::binomial())
  fit_c <- stats::glm(c_flag ~ A   + L_1 + k,
                      data = pt, family = stats::binomial())

  cy <- stats::coef(fit_y); cd <- stats::coef(fit_d); cc <- stats::coef(fit_c)
  # Use absolute tolerance — relative tolerance from expect_equal punishes
  # small true coefficients. n=8000 + rare events (especially censoring)
  # gives Wald SE on the order of 0.05–0.15 per coefficient.
  abs_tol <- 0.20

  expect_lt(abs(unname(cy["(Intercept)"]) - beta_y_true$intercept), abs_tol)
  expect_lt(abs(unname(cy["A_y"])         - beta_y_true$A_Y),       abs_tol)
  expect_lt(abs(unname(cy["L_1"])         - beta_y_true$L),         abs_tol)
  expect_lt(abs(unname(cy["k"])           - beta_y_true$k),         abs_tol)

  expect_lt(abs(unname(cd["(Intercept)"]) - beta_d_true$intercept), abs_tol)
  expect_lt(abs(unname(cd["A_d"])         - beta_d_true$A_D),       abs_tol)
  expect_lt(abs(unname(cd["L_1"])         - beta_d_true$L),         abs_tol)

  expect_lt(abs(unname(cc["(Intercept)"]) - beta_c_true$intercept), abs_tol)
  expect_lt(abs(unname(cc["A"])           - beta_c_true$A),         abs_tol)
  expect_lt(abs(unname(cc["L_1"])         - beta_c_true$L),         abs_tol)
})


# ============================================================================
# Step 2: discrete S(k) matches KM with NO censoring
# ============================================================================
# Y-only DGP — h_D and h_C are effectively zero (intercept = -50). Only
# event types {0, 1} appear: Y events (1) and admin-censoring at K (0).
# Discrete S(t) from the pooled logistic should match KM at all integer t.
# ============================================================================

test_that("Step 2: discrete S(k) from pooled logistic matches KM (no censoring)", {
  skip_on_cran()
  if (!requireNamespace("survival", quietly = TRUE)) skip("survival not available")

  sim <- simulate_separable_data(
    n      = 5000, K = 12, prev_L = numeric(0),
    beta_y = list(intercept = -2.5, L = numeric(0), A_Y = 0, A_D = 0, k = 0),
    beta_d = list(intercept = -50,  L = numeric(0), A_Y = 0, A_D = 0, k = 0),
    beta_c = list(intercept = -50,  L = numeric(0), A   = 0,         k = 0),
    beta_a = list(intercept = 0,    L = numeric(0)),
    seed   = 2
  )
  expect_setequal(unique(sim$data$event_type), c(0L, 1L))

  pt <- build_pt(sim$data, K = 12)
  fit_y <- stats::glm(y_flag ~ k + I(k^2), data = pt, family = stats::binomial())

  k_seq   <- 0:11
  h_y_hat <- stats::predict(fit_y, newdata = data.frame(k = k_seq),
                            type = "response")
  # S_pkg(t) = prod_{j<t}(1 - h_y(j)); evaluated at t = 0, 1, ..., 12
  S_pkg <- c(1, cumprod(1 - h_y_hat))

  km <- survival::survfit(
    survival::Surv(time = event_time + 1L,
                   event = (event_type == 1L)) ~ 1,
    data = sim$data
  )
  km_summary <- summary(km, times = 0:12, extend = TRUE)
  S_km <- km_summary$surv

  expect_equal(length(S_pkg), length(S_km))
  expect_lt(max(abs(S_pkg - S_km)), 0.03)
})


# ============================================================================
# Step 3: discrete S(k) matches KM with INDEPENDENT censoring
# ============================================================================
# Same as step 2 but adds independent censoring (h_C does not depend on A
# or L). KM is unbiased for S_Y when censoring is independent, and so is
# our pooled-logistic + cumulative-product estimator.
# ============================================================================

test_that("Step 3: discrete S(k) from pooled logistic matches KM (indep censoring)", {
  skip_on_cran()
  if (!requireNamespace("survival", quietly = TRUE)) skip("survival not available")

  sim <- simulate_separable_data(
    n      = 5000, K = 12, prev_L = numeric(0),
    beta_y = list(intercept = -2.5, L = numeric(0), A_Y = 0, A_D = 0, k = 0),
    beta_d = list(intercept = -50,  L = numeric(0), A_Y = 0, A_D = 0, k = 0),
    beta_c = list(intercept = -3.5, L = numeric(0), A   = 0,         k = 0),
    beta_a = list(intercept = 0,    L = numeric(0)),
    seed   = 3
  )
  expect_true(all(sim$data$event_type %in% c(0L, 1L)))

  pt <- build_pt(sim$data, K = 12)
  fit_y <- stats::glm(y_flag ~ k + I(k^2), data = pt, family = stats::binomial())

  k_seq   <- 0:11
  h_y_hat <- stats::predict(fit_y, newdata = data.frame(k = k_seq),
                            type = "response")
  S_pkg <- c(1, cumprod(1 - h_y_hat))

  km <- survival::survfit(
    survival::Surv(time = event_time + 1L,
                   event = (event_type == 1L)) ~ 1,
    data = sim$data
  )
  km_summary <- summary(km, times = 0:12, extend = TRUE)
  S_km <- km_summary$surv

  expect_equal(length(S_pkg), length(S_km))
  expect_lt(max(abs(S_pkg - S_km)), 0.03)
})
