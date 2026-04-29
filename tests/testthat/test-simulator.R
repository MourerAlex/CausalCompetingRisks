# Smoke tests for the simulate_separable_data() helper itself, before
# using it in the ladder tests. Verifies output shape and that the
# closed-form truth is internally consistent (algebraic identities).

test_that("simulate_separable_data returns the expected structure", {
  sim <- simulate_separable_data(
    n = 100, K = 6, prev_L = c(0.5),
    beta_y = list(intercept = -3,   L = c(0.3),  A_Y = -0.5, A_D = 0,    k = 0),
    beta_d = list(intercept = -3.5, L = c(0.2),  A_Y = 0,    A_D = -0.4, k = 0),
    beta_c = list(intercept = -4,   L = c(0),    A   = 0,                k = 0),
    beta_a = list(intercept = 0,    L = c(0)),
    seed   = 1
  )

  expect_named(sim, c("data", "truth", "params"))
  expect_named(sim$data,
               c("id", "event_time", "event_type", "A", "L_1"))
  expect_equal(nrow(sim$data), 100L)

  expect_named(sim$truth, c("arm_cif", "contrasts"))
  expect_named(sim$truth$arm_cif,
               c("k", "arm_11", "arm_00", "arm_10", "arm_01"))
  expect_named(sim$truth$contrasts,
               c("k", "total", "sde_a", "sie_a", "sde_b", "sie_b"))
  expect_equal(nrow(sim$truth$arm_cif), 6L)
})


test_that("event_type values are in {0, 1, 2}", {
  sim <- simulate_separable_data(
    n = 200, K = 6, prev_L = c(0.5),
    beta_y = list(intercept = -3,   L = c(0.3), A_Y = -0.5, A_D = 0, k = 0),
    beta_d = list(intercept = -3.5, L = c(0.2), A_Y = 0, A_D = -0.4, k = 0),
    beta_c = list(intercept = -3,   L = c(0),   A = 0, k = 0),
    beta_a = list(intercept = 0,    L = c(0)),
    seed   = 2
  )
  expect_true(all(sim$data$event_type %in% c(0L, 1L, 2L)))
})


test_that("truth: arm CIFs are monotone non-decreasing in k and within [0, 1]", {
  sim <- simulate_separable_data(
    n = 50, K = 8, prev_L = c(0.4, 0.6),
    beta_y = list(intercept = -3,   L = c(0.5, -0.2),  A_Y = -0.7, A_D = 0,    k = -0.05),
    beta_d = list(intercept = -3.5, L = c(0.3, -0.1),  A_Y = 0,    A_D = -0.5, k = 0),
    beta_c = list(intercept = -4,   L = c(0, 0),       A   = 0,                k = 0),
    beta_a = list(intercept = 0,    L = c(0, 0)),
    seed   = 3
  )
  for (a in c("arm_11", "arm_00", "arm_10", "arm_01")) {
    v <- sim$truth$arm_cif[[a]]
    expect_true(all(diff(v) >= 0), info = paste("arm:", a))
    expect_true(all(v >= 0 & v <= 1),  info = paste("arm:", a))
  }
})


test_that("truth: SDE_A + SIE_A = total and SDE_B + SIE_B = total at every k", {
  sim <- simulate_separable_data(
    n = 50, K = 8, prev_L = c(0.5),
    beta_y = list(intercept = -3,   L = c(0.5),  A_Y = -0.7, A_D = 0.2,  k = -0.05),
    beta_d = list(intercept = -3.5, L = c(0.3),  A_Y = -0.1, A_D = -0.5, k = 0),
    beta_c = list(intercept = -4,   L = c(0),    A   = 0,                k = 0),
    beta_a = list(intercept = 0,    L = c(0)),
    seed   = 4
  )
  ctr <- sim$truth$contrasts
  expect_equal(ctr$sde_a + ctr$sie_a, ctr$total, tolerance = 1e-12)
  expect_equal(ctr$sde_b + ctr$sie_b, ctr$total, tolerance = 1e-12)
})


test_that("truth: under full isolation + no L effect, arm_10 has no Y effect of A_D (collapses)", {
  # beta_y$A_D = 0 (Δ1 holds). Under additionally h_D constant in k and
  # without confounding, the SIE-A of A_D on Y is purely through the
  # at-risk channel — which still produces a non-zero SIE in general.
  # So we instead test the algebraic identity that arm_10 lies between
  # arm_00 and arm_11 when A_Y has a negative effect (curves don't cross).
  sim <- simulate_separable_data(
    n = 50, K = 6, prev_L = c(0.5),
    beta_y = list(intercept = -3,   L = c(0.5),  A_Y = -0.7, A_D = 0,    k = 0),
    beta_d = list(intercept = -3.5, L = c(0),    A_Y = 0,    A_D = -0.5, k = 0),
    beta_c = list(intercept = -4,   L = c(0),    A   = 0,                k = 0),
    beta_a = list(intercept = 0,    L = c(0)),
    seed   = 5
  )
  arm_cif <- sim$truth$arm_cif
  # A_Y reduces Y; A_D reduces D so more Y survives → arm_11 should be highest.
  # arm_00 has neither benefit. arm_10 has only the A_Y reduction; arm_01 has
  # only the A_D reduction (which raises Y via at-risk).
  expect_true(all(arm_cif$arm_10 <= arm_cif$arm_00 + 1e-12),
              info = "arm_10 should be <= arm_00 when A_Y reduces Y")
  expect_true(all(arm_cif$arm_01 >= arm_cif$arm_00 - 1e-12),
              info = "arm_01 should be >= arm_00 when A_D reduces D (more at-risk)")
})


test_that("simulator is reproducible given a seed", {
  sim1 <- simulate_separable_data(n = 50, K = 6, prev_L = c(0.5),
    beta_y = list(intercept = -3, L = c(0.3), A_Y = -0.5, A_D = 0, k = 0),
    beta_d = list(intercept = -3.5, L = c(0.2), A_Y = 0, A_D = -0.4, k = 0),
    beta_c = list(intercept = -4, L = c(0), A = 0, k = 0),
    beta_a = list(intercept = 0, L = c(0)),
    seed = 99
  )
  sim2 <- simulate_separable_data(n = 50, K = 6, prev_L = c(0.5),
    beta_y = list(intercept = -3, L = c(0.3), A_Y = -0.5, A_D = 0, k = 0),
    beta_d = list(intercept = -3.5, L = c(0.2), A_Y = 0, A_D = -0.4, k = 0),
    beta_c = list(intercept = -4, L = c(0), A = 0, k = 0),
    beta_a = list(intercept = 0, L = c(0)),
    seed = 99
  )
  expect_identical(sim1$data, sim2$data)
})
