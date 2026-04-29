# Tests for assumptions() accessor and S3 print/format methods
#
# Note: the `$isolation` slot (numeric swap-weight read-out) is currently
# omitted from assumptions() pending a math check against the Stensrud
# Appendix. Tests for that slot have been removed; restore alongside the
# slot when the math link is confirmed.

# ---- Fixture: fit with IPW ----
make_fit_ipw <- function() {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("normal_act"),
    event_y = 1, event_d = 2, event_c = 0, n_intervals = 4
  )
  suppressWarnings(separable_effects(pt, method = c("ipw1")))
}

make_fit_gformula_only <- function() {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("normal_act"),
    event_y = 1, event_d = 2, event_c = 0, n_intervals = 4
  )
  suppressWarnings(separable_effects(pt, method = c("gformula")))
}


test_that("assumptions() returns object of class separable_effects_assumptions", {
  fit <- make_fit_ipw()
  a <- assumptions(fit)
  expect_s3_class(a, "separable_effects_assumptions")
})


test_that("assumptions() has the locked slot structure (isolation deferred)", {
  fit <- make_fit_ipw()
  a <- assumptions(fit)
  expect_named(
    a,
    c("estimand", "assumptions",
      "dismissible_D1", "dismissible_D2")
  )
})


test_that("estimand slot has target, description, arms, two decompositions", {
  fit <- make_fit_ipw()
  est <- assumptions(fit)$estimand
  expect_named(est,
               c("target", "description", "arms",
                 "decomposition_A", "decomposition_B"))
  expect_length(est$arms, 4L)
  expect_named(est$decomposition_A, c("label", "sde", "sie"))
  expect_named(est$decomposition_B, c("label", "sie", "sde"))
})


test_that("each assumption entry has statement, testable, citation, diagnostic", {
  fit <- make_fit_ipw()
  a <- assumptions(fit)$assumptions
  for (nm in names(a)) {
    expect_true(all(c("statement", "testable", "citation", "diagnostic")
                    %in% names(a[[nm]])),
                info = paste("entry:", nm))
  }
})


test_that("standard assumptions list includes GDA, E1, E2, consistency, positivity", {
  fit <- make_fit_ipw()
  a <- assumptions(fit)$assumptions
  expect_true(all(c("GDA",
                    "E1_treatment_exchangeability",
                    "E2_censoring_exchangeability",
                    "consistency",
                    "positivity") %in% names(a)))
})


test_that("dismissible_D1 and dismissible_D2 share the same field shape", {
  fit <- make_fit_ipw()
  a <- assumptions(fit)
  for (nm in c("dismissible_D1", "dismissible_D2")) {
    entry <- a[[nm]]
    expect_true(all(c("statement", "testable", "citation", "diagnostic")
                    %in% names(entry)),
                info = paste("entry:", nm))
  }
})


# Isolation slot tests removed pending math check (see header note).


test_that("print.separable_effects_assumptions runs without error", {
  fit <- make_fit_ipw()
  a <- assumptions(fit)
  expect_output(print(a), "Identifying Assumptions")
  expect_output(print(a), "Estimand")
  expect_output(print(a), "GDA")
  expect_output(print(a), "dismissible_D1")
  expect_output(print(a), "dismissible_D2")
})


test_that("format(x, style = 'markdown') returns a single string with headers", {
  fit <- make_fit_ipw()
  md <- format(assumptions(fit), style = "markdown")
  expect_type(md, "character")
  expect_length(md, 1L)
  expect_match(md, "## Estimand", fixed = TRUE)
  expect_match(md, "## Identifying assumptions", fixed = TRUE)
})


test_that("format() rejects unknown styles", {
  fit <- make_fit_ipw()
  expect_error(
    format(assumptions(fit), style = "html"),
    "Only style = \"markdown\""
  )
})


test_that("print.separable_effects shows one-line identification header", {
  fit <- make_fit_ipw()
  expect_output(print(fit), "Estimand")
  expect_output(print(fit), "assumptions\\(fit\\)")
})
