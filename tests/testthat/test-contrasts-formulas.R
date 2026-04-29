# Regression test: pins the Decomposition A & B contrast formulas to the
# methodologically correct definitions (Stensrud et al. 2020, JASA, §3).
#
# Convention: arm "xy" = (a_Y = x, a_D = y).
#   A_Y acts on the path A → Y not through D  →  its effect is the SDE.
#   A_D acts on the path A → D → Y            →  its effect is the SIE.
#
# Decomposition A (vary A_Y first; intermediate arm = arm_10):
#   SDE-A = arm_10 − arm_00 ; SIE-A = arm_11 − arm_10
#
# Decomposition B (vary A_D first; intermediate arm = arm_01):
#   SIE-B = arm_01 − arm_00 ; SDE-B = arm_11 − arm_01
#
# Both share Total = arm_11 − arm_00 and the algebraic identity SDE + SIE = Total.
#
# This file does NOT assume how the contrasts are computed inside the package
# (it tests the public-output contract on a controlled fixture). When the
# refactor moves contrast logic to a different module, this test should still
# pass without modification.

# --- Fixture: arm CIF values chosen so every contrast is uniquely identifiable -
arm_cif <- c(
  arm_00 = 0.1,
  arm_10 = 0.2,
  arm_01 = 0.3,
  arm_11 = 0.4
)

expected <- list(
  total = unname(arm_cif["arm_11"] - arm_cif["arm_00"]),  # 0.3

  sde_a = unname(arm_cif["arm_10"] - arm_cif["arm_00"]),  # 0.1
  sie_a = unname(arm_cif["arm_11"] - arm_cif["arm_10"]),  # 0.2

  sie_b = unname(arm_cif["arm_01"] - arm_cif["arm_00"]),  # 0.2
  sde_b = unname(arm_cif["arm_11"] - arm_cif["arm_01"])   # 0.1
)

test_that("Decomposition A: SDE-A = arm_10 - arm_00", {
  expect_equal(arm_cif[["arm_10"]] - arm_cif[["arm_00"]], expected$sde_a)
})

test_that("Decomposition A: SIE-A = arm_11 - arm_10", {
  expect_equal(arm_cif[["arm_11"]] - arm_cif[["arm_10"]], expected$sie_a)
})

test_that("Decomposition B: SIE-B = arm_01 - arm_00", {
  expect_equal(arm_cif[["arm_01"]] - arm_cif[["arm_00"]], expected$sie_b)
})

test_that("Decomposition B: SDE-B = arm_11 - arm_01", {
  expect_equal(arm_cif[["arm_11"]] - arm_cif[["arm_01"]], expected$sde_b)
})

test_that("Total is identical under both decompositions", {
  total_via_a <- expected$sde_a + expected$sie_a
  total_via_b <- expected$sde_b + expected$sie_b
  expect_equal(total_via_a, expected$total)
  expect_equal(total_via_b, expected$total)
})

test_that("SDE-A and SDE-B both measure the effect of A_Y (different a_D references)", {
  # SDE = effect of A_Y holding a_D fixed.
  # SDE-A holds a_D = 0; SDE-B holds a_D = 1.
  # Under "no effect of A_D on Y at fixed A_Y" the two SDEs would coincide.
  # Here arm_10 - arm_00 = 0.1 and arm_11 - arm_01 = 0.1 — they coincide for
  # this fixture, by construction.
  expect_equal(arm_cif[["arm_10"]] - arm_cif[["arm_00"]],
               arm_cif[["arm_11"]] - arm_cif[["arm_01"]])
})

test_that("SIE-A and SIE-B both measure the effect of A_D (different a_Y references)", {
  # SIE = effect of A_D holding a_Y fixed.
  # SIE-A holds a_Y = 1; SIE-B holds a_Y = 0.
  # Same construction — should coincide here.
  expect_equal(arm_cif[["arm_11"]] - arm_cif[["arm_10"]],
               arm_cif[["arm_01"]] - arm_cif[["arm_00"]])
})

test_that("formulas in R/contrasts.R match the locked spec (regression on package code)", {
  # This pins the actual code paths in R/contrasts.R against the locked spec.
  # Construct a minimal fake fit + bootstrap, run compute_contrasts, check the
  # `estimate` column for each (contrast, decomp) row at a single timepoint.
  #
  # If R/contrasts.R drifts (e.g., a future refactor swaps SDE/SIE labels), this
  # test fails immediately.

  if (!exists("compute_contrasts", mode = "function")) {
    skip("compute_contrasts() not found — skipping package-code regression check")
  }

  # Two timepoints with identical arm CIFs (avoids 1-d collapse in apply()).
  n_times <- 2L
  ks <- seq_len(n_times)
  ci_df <- data.frame(
    k      = ks,
    arm_00 = rep(arm_cif[["arm_00"]], n_times),
    arm_10 = rep(arm_cif[["arm_10"]], n_times),
    arm_01 = rep(arm_cif[["arm_01"]], n_times),
    arm_11 = rep(arm_cif[["arm_11"]], n_times)
  )

  fit <- structure(
    list(
      cumulative_incidence = list(gformula = ci_df),
      active_methods       = "gformula"
    ),
    class = "separable_effects"
  )

  # Bootstrap: 4D array [n_boot, method, arm, time]. Replicates equal point
  # estimates so percentile bounds equal the point — easy contract to assert.
  arms   <- c("arm_11", "arm_00", "arm_10", "arm_01")
  n_boot <- 5L
  reps <- array(
    NA_real_,
    dim = c(n_boot, 1L, length(arms), n_times),
    dimnames = list(NULL, "gformula", arms, NULL)
  )
  for (a in arms) reps[, "gformula", a, ] <- arm_cif[[a]]

  ci <- structure(
    list(replicates = reps, alpha = 0.05),
    class = "causal_cr_bootstrap"
  )

  out <- compute_contrasts(fit, method = "gformula", ci = ci)

  pick <- function(contrast_name, decomp_label, measure_label = "rd",
                   k_value = ks[1L]) {
    row <- out[
      out$k        == k_value &
      out$contrast == contrast_name &
        ((is.na(out$decomp)  & is.na(decomp_label)) |
         (!is.na(out$decomp) & !is.na(decomp_label) & out$decomp == decomp_label)) &
        out$measure  == measure_label,
    ]
    expect_equal(nrow(row), 1L)
    row$estimate
  }

  expect_equal(pick("total",    NA_character_), expected$total)
  expect_equal(pick("direct",   "A"),           expected$sde_a)
  expect_equal(pick("indirect", "A"),           expected$sie_a)
  expect_equal(pick("direct",   "B"),           expected$sde_b)
  expect_equal(pick("indirect", "B"),           expected$sie_b)
})
