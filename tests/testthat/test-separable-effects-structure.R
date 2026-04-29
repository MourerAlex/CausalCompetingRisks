# Structural tests for separable_effects() return object

# ---- Helper: minimal realistic pt fixture ----
make_pt <- function(n = 30, seed = 1) {
  set.seed(seed)
  df <- data.frame(
    id         = seq_len(n),
    event_time = sample(3:30, n, replace = TRUE),
    event_type = sample(c(0L, 1L, 2L), n, replace = TRUE,
                        prob = c(0.3, 0.35, 0.35)),
    A          = sample(c(0L, 1L), n, replace = TRUE),
    cov1       = sample(c(0L, 1L), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  suppressWarnings(to_person_time(
    df,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0,
    n_intervals = 6
  ))
}


# ---- Return structure ----

test_that("separable_effects returns S3 object of class 'separable_effects'", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt))
  expect_s3_class(fit, "separable_effects")
})


test_that("method = 'all' produces all three cumulative_incidence entries", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "all"))
  expect_true(is.list(fit$cumulative_incidence))
  expect_setequal(names(fit$cumulative_incidence),
                  c("gformula", "ipw1", "ipw2"))
})


test_that("method = 'gformula' produces only gformula entry", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "gformula"))
  expect_equal(names(fit$cumulative_incidence), "gformula")
})


test_that("method = 'ipw1' produces only ipw1 entry", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "ipw1"))
  expect_equal(names(fit$cumulative_incidence), "ipw1")
})


test_that("method = 'ipw2' produces only ipw2 entry", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "ipw2"))
  expect_equal(names(fit$cumulative_incidence), "ipw2")
})


# ---- Arm columns ----

test_that("each cumulative_incidence data.frame has k and all four arms", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "all"))
  for (m in names(fit$cumulative_incidence)) {
    df <- fit$cumulative_incidence[[m]]
    # All methods (g-formula, ipw1, ipw2) emit all 4 arms including arm_01
    # (Decomposition B sensitivity).
    expect_true(all(c("k", "arm_11", "arm_00", "arm_10", "arm_01") %in% names(df)),
                info = paste("method:", m))
    # Cumulative incidence should be non-decreasing for every arm
    for (a in c("arm_11", "arm_00", "arm_10", "arm_01")) {
      expect_true(all(diff(df[[a]]) >= -1e-9),
                  info = paste("method:", m, "arm:", a, "not monotone"))
    }
  }
})


# ---- weights slot ----

test_that("gformula-only method leaves weights NULL", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "gformula"))
  expect_null(fit$weights)
})


test_that("IPW methods populate weights with expected components", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "all"))
  expect_false(is.null(fit$weights))
  expect_true(all(c("pt_data_weighted", "weight_summary", "flagged_ids",
                    "truncate")
                  %in% names(fit$weights)))
})


test_that("pt_data_weighted carries raw and truncated weight columns", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "all"))
  w <- fit$weights$pt_data_weighted
  expect_true("w_cens_raw" %in% names(w))
  expect_true("w_cens" %in% names(w))
  expect_true("w_d_arm_10_raw" %in% names(w))
  expect_true("w_d_arm_10" %in% names(w))
  expect_true("w_d_arm_01_raw" %in% names(w))
  expect_true("w_d_arm_01" %in% names(w))
})


test_that("pt_data_weighted carries hazard and cumprod columns", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "all"))
  w <- fit$weights$pt_data_weighted
  expected_cols <- c(
    "haz_y_a1", "haz_y_a0", "haz_d_a1", "haz_d_a0", "haz_c",
    "cumprod_one_minus_hazd_a1", "cumprod_one_minus_hazd_a0",
    "cumprod_one_minus_hazy_a1", "cumprod_one_minus_hazy_a0",
    "lag_cumprod_one_minus_hazy_a1", "lag_cumprod_one_minus_hazy_a0"
  )
  for (col in expected_cols) {
    expect_true(col %in% names(w),
                info = paste("missing:", col))
  }
})


# ---- ipcw = FALSE ----

test_that("ipcw = FALSE skips model_c", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "all", ipcw = FALSE))
  expect_null(fit$models$model_c)
})


test_that("ipcw = FALSE omits w_cens columns", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "all", ipcw = FALSE))
  if (!is.null(fit$weights)) {
    w <- fit$weights$pt_data_weighted
    expect_false("w_cens" %in% names(w))
    expect_false("haz_c" %in% names(w))
  }
})


# ---- warnings slot ----

test_that("warnings are captured into fit$warnings", {
  # Force a constant covariate to trigger a validation warning that
  # propagates through the withCallingHandlers wrapper.
  df <- data.frame(
    id         = 1:30,
    event_time = sample(3:30, 30, replace = TRUE),
    event_type = sample(c(0L, 1L, 2L), 30, replace = TRUE),
    A          = sample(c(0L, 1L), 30, replace = TRUE),
    cov1       = 1L,  # constant
    stringsAsFactors = FALSE
  )
  pt <- suppressWarnings(to_person_time(
    df,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0,
    n_intervals = 6
  ))
  fit <- suppressWarnings(separable_effects(pt, method = "gformula"))
  expect_type(fit$warnings, "character")
})


# ---- method validation ----

test_that("unknown method errors early", {
  pt <- make_pt()
  expect_error(
    suppressWarnings(separable_effects(pt, method = "nonsense")),
    "Unknown method"
  )
})


# ---- n, times, models slots ----

test_that("fit exposes n, times, models, call", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "gformula"))
  expect_equal(fit$n, 30)
  expect_equal(fit$times, attr(pt, "cut_times"))
  expect_true(!is.null(fit$models$model_y))
  expect_true(!is.null(fit$models$model_d))
  expect_true(!is.null(fit$call))
})


# ---- method as character vector ----

test_that("method accepts a character vector subset", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = c("gformula", "ipw2")))
  expect_setequal(names(fit$cumulative_incidence),
                  c("gformula", "ipw2"))
})


test_that("method with unknown entry errors", {
  pt <- make_pt()
  expect_error(
    suppressWarnings(separable_effects(pt, method = c("gformula", "nonsense"))),
    "Unknown method"
  )
})


# ---- model_checks slot ----

test_that("model_checks slot is a list with expected structure", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "all"))
  expect_true(is.list(fit$model_checks))
  # Should have entries for any fitted model
  expect_true(!is.null(fit$model_checks$y))
  expect_true(!is.null(fit$model_checks$d))
})


test_that("model_checks$y has the expected fields", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "gformula"))
  chk <- fit$model_checks$y
  expect_equal(chk$label, "Y-hazard")
  expect_type(chk$converged, "logical")
  expect_type(chk$min_fitted, "double")
  expect_type(chk$max_fitted, "double")
  expect_type(chk$glm_warnings, "character")
})


test_that("gformula-only leaves model_checks$c NULL", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "gformula"))
  expect_null(fit$model_checks$c)
})


test_that("ipcw = FALSE leaves model_checks$c NULL", {
  pt <- make_pt()
  fit <- suppressWarnings(separable_effects(pt, method = "all",
                                      ipcw = FALSE))
  expect_null(fit$model_checks$c)
})
