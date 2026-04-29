test_that("fit_separable_effects returns plain list (no class)", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A",
    covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
    event_y = 1, event_d = 2, event_c = 0
  )

  out <- suppressWarnings(fit_separable_effects(
    pt_data = pt,
    id_col = "id",
    treatment_col = "A",
    covariates_vec = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
    cut_times = attr(pt, "cut_times"),
    active_methods = c("gformula"),
    formulas = NULL,
    ipcw = TRUE,
    truncate = c(0.01, 0.99)
  ))

  expect_false(inherits(out, "separable_effects"))
  expect_true(is.list(out))
  expect_named(out, c("cumulative_incidence", "weights", "models", "model_checks"))
  expect_true("gformula" %in% names(out$cumulative_incidence))
})


test_that("separable_effects stores bootstrap-needed args on fit", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A",
    covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
    event_y = 1, event_d = 2, event_c = 0
  )

  fit <- suppressWarnings(separable_effects(pt, method = c("gformula")))

  expect_equal(fit$treatment_col, "A")
  expect_equal(fit$active_methods, "gformula")
  expect_equal(fit$covariates, c("normal_act", "age_cat", "cv_hist", "hemo_bin"))
  expect_equal(fit$ipcw, TRUE)
  expect_equal(fit$truncate, c(0.01, 0.99))
})


test_that("bootstrap returns separable_effects_bootstrap with 4D replicates", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A",
    covariates = c("normal_act"),   # keep tiny for speed
    event_y = 1, event_d = 2, event_c = 0,
    n_intervals = 4
  )

  fit <- suppressWarnings(separable_effects(pt, method = c("gformula")))

  # Tiny bootstrap for speed
  boot <- suppressMessages(suppressWarnings(bootstrap(fit, n_boot = 5)))

  expect_s3_class(boot, "separable_effects_bootstrap")
  expect_equal(length(dim(boot$replicates)), 4)
  expect_equal(dim(boot$replicates)[1], 5)         # n_boot
  expect_equal(dim(boot$replicates)[2], 1)         # n_methods (just gformula)
  expect_equal(dim(boot$replicates)[3], 4)         # n_arms (includes arm_01)
  expect_equal(dimnames(boot$replicates)[[2]], "gformula")
  expect_equal(dimnames(boot$replicates)[[3]],
               c("arm_11", "arm_00", "arm_10", "arm_01"))
})


test_that("bootstrap ci_curves has lower/upper per arm per method", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("normal_act"),
    event_y = 1, event_d = 2, event_c = 0, n_intervals = 4
  )

  fit <- suppressWarnings(separable_effects(pt, method = c("gformula")))
  boot <- suppressMessages(suppressWarnings(bootstrap(fit, n_boot = 5)))

  expect_named(boot$ci_curves, "gformula")
  gci <- boot$ci_curves$gformula
  expect_true(all(c("k",
                    "arm_11_lower", "arm_11_upper",
                    "arm_00_lower", "arm_00_upper",
                    "arm_10_lower", "arm_10_upper",
                    "arm_01_lower", "arm_01_upper") %in% names(gci)))
})


test_that("contrast() errors without ci argument", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("normal_act"),
    event_y = 1, event_d = 2, event_c = 0, n_intervals = 4
  )
  fit <- suppressWarnings(separable_effects(pt, method = c("gformula")))

  expect_error(contrast(fit, method = "gformula"), "'ci' argument is required")
  expect_error(contrast(fit), "'ci' argument is required")
})


test_that("contrast() defaults to max time (10 rows: 1 total + 2 direct + 2 indirect, RD/RR)", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("normal_act"),
    event_y = 1, event_d = 2, event_c = 0, n_intervals = 4
  )
  fit <- suppressWarnings(separable_effects(pt, method = c("gformula")))
  boot <- suppressMessages(suppressWarnings(bootstrap(fit, n_boot = 5)))

  ctr <- contrast(fit, method = "gformula", ci = boot)

  expect_s3_class(ctr, "separable_effects_contrast")
  expect_true(is.data.frame(ctr$contrasts))
  expect_named(ctr$contrasts,
               c("k", "contrast", "decomp", "measure",
                 "estimate", "lower", "upper"))
  # 1 time point x 10 contrast rows
  expect_equal(nrow(ctr$contrasts), 10)
  # Default time = max
  expect_equal(unique(ctr$contrasts$k), max(fit$times))
  expect_equal(ctr$time, max(fit$times))
  # Decomp values
  expect_equal(sort(unique(ctr$contrasts$contrast)),
               c("sep_direct", "sep_indirect", "total"))
  expect_true(all(ctr$contrasts$decomp[ctr$contrasts$contrast == "total"] |>
                    is.na()))
  expect_true(all(c("A", "B") %in% ctr$contrasts$decomp))
  expect_equal(sort(unique(ctr$contrasts$measure)), c("rd", "rr"))
})


test_that("contrast() with time = <value> snaps to nearest cut and filters", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("normal_act"),
    event_y = 1, event_d = 2, event_c = 0, n_intervals = 4
  )
  fit <- suppressWarnings(separable_effects(pt, method = c("gformula")))
  boot <- suppressMessages(suppressWarnings(bootstrap(fit, n_boot = 5)))

  # Pick the first cut time, no snapping
  k_target <- fit$times[1]
  ctr <- contrast(fit, method = "gformula", ci = boot, time = k_target)
  expect_equal(nrow(ctr$contrasts), 10)
  expect_equal(unique(ctr$contrasts$k), k_target)
  expect_equal(ctr$time, k_target)

  # An off-cut value should snap and emit a message
  off_cut <- mean(fit$times[1:2])
  expect_message(
    contrast(fit, method = "gformula", ci = boot, time = off_cut),
    "snapped to nearest cut time"
  )
})


test_that("contrast() for IPW has populated decomp B (IPW now emits arm_01)", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("normal_act"),
    event_y = 1, event_d = 2, event_c = 0, n_intervals = 4
  )
  fit <- suppressWarnings(separable_effects(pt, method = c("ipw")))
  boot <- suppressMessages(suppressWarnings(bootstrap(fit, n_boot = 5)))

  ctr <- contrast(fit, method = "ipw_rep1", ci = boot)
  # Decomp B rows should now be populated (non-NA) because IPW Rep 1
  # dispatcher emits arm_01 alongside the other three arms.
  decomp_B <- ctr$contrasts[!is.na(ctr$contrasts$decomp) &
                             ctr$contrasts$decomp == "B", ]
  expect_true(nrow(decomp_B) > 0)
  expect_true(all(!is.na(decomp_B$estimate)))
})


test_that("risk() returns long-format $risk; lower/upper populated only with ci", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("normal_act"),
    event_y = 1, event_d = 2, event_c = 0, n_intervals = 4
  )
  fit <- suppressWarnings(separable_effects(pt, method = c("gformula")))
  boot <- suppressMessages(suppressWarnings(bootstrap(fit, n_boot = 5)))

  r_no_ci <- risk(fit)
  expect_true(is.data.frame(r_no_ci$risk))
  expect_named(r_no_ci$risk,
               c("method", "arm", "a_y", "a_d", "k", "value", "lower", "upper"))
  expect_true(all(is.na(r_no_ci$risk$lower)))
  expect_true(all(is.na(r_no_ci$risk$upper)))

  r_with_ci <- risk(fit, ci = boot)
  expect_true(is.data.frame(r_with_ci$risk))
  expect_false(all(is.na(r_with_ci$risk$lower)))
  expect_false(all(is.na(r_with_ci$risk$upper)))
})


test_that("separable_effects errors on NULL / non-data.frame / empty input", {
  expect_error(separable_effects(NULL), "pt_data is NULL")
  expect_error(separable_effects(matrix(1:4, 2, 2)), "must be a data.frame")
  expect_error(separable_effects(list(a = 1)), "must be a data.frame")
  expect_error(
    separable_effects(data.frame(id = integer(0), k = integer(0),
                         A = integer(0), y_flag = integer(0),
                         d_flag = integer(0), c_flag = integer(0)),
              id = "id", treatment = "A"),
    "0 rows"
  )
})


test_that("to_person_time errors on NULL / non-data.frame / empty input", {
  expect_error(
    to_person_time(NULL, id = "id", time = "event_time", event = "event_type",
                   treatment = "A", covariates = character(),
                   event_y = 1, event_d = 2, event_c = 0),
    "data is NULL"
  )
  expect_error(
    to_person_time(matrix(1:4, 2, 2), id = "id", time = "event_time",
                   event = "event_type", treatment = "A",
                   covariates = character(),
                   event_y = 1, event_d = 2, event_c = 0),
    "must be a data.frame"
  )
  expect_error(
    to_person_time(
      data.frame(id = integer(0), event_time = integer(0),
                 event_type = integer(0), A = integer(0)),
      id = "id", time = "event_time", event = "event_type",
      treatment = "A", covariates = character(),
      event_y = 1, event_d = 2, event_c = 0
    ),
    "0 rows"
  )
})


test_that("separable_effects errors if id/treatment passed on a person_time", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("normal_act"),
    event_y = 1, event_d = 2, event_c = 0, n_intervals = 4
  )

  expect_error(
    suppressWarnings(separable_effects(pt, id = "id")),
    "Do not pass them again"
  )
  expect_error(
    suppressWarnings(separable_effects(pt, treatment = "A")),
    "Do not pass them again"
  )
  expect_error(
    suppressWarnings(separable_effects(pt, covariates = c("normal_act"))),
    "Do not pass them again"
  )
})


test_that("bootstrap validates n_boot and alpha", {
  pt <- to_person_time(
    prostate_data,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("normal_act"),
    event_y = 1, event_d = 2, event_c = 0, n_intervals = 4
  )
  fit <- suppressWarnings(separable_effects(pt, method = c("gformula")))

  expect_error(bootstrap(fit, n_boot = 0), "positive integer")
  expect_error(bootstrap(fit, n_boot = 1.5), "positive integer")
  expect_error(bootstrap(fit, n_boot = 5, alpha = 0), "in \\(0, 1\\)")
  expect_error(bootstrap(fit, n_boot = 5, alpha = 1.5), "in \\(0, 1\\)")
})
