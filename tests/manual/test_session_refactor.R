# =============================================================================
# Manual test of everything reviewed this session
# =============================================================================
# Source files and load data, then step through each section and inspect
# results. Each block prints a header and a tidy summary.

# ---- Setup ----
setwd("/home/moureralex/Bureau/cowork/separable_effects")

source("R/utils.R")
source("R/validate.R")
source("R/data_prep.R")
source("R/hazard_models.R")
source("R/gformula.R")
source("R/ipw.R")
source("R/contrasts.R")
source("R/bootstrap.R")
source("R/accessors.R")
source("R/print.R")
source("R/plot.R")
source("R/separable_effects.R")

load("data/prostate_data.rda")


# =============================================================================
# 1. check_covariate_quality refactor
# =============================================================================
cat("\n============================================================\n")
cat(" 1. check_covariate_quality\n")
cat("============================================================\n\n")

cat("1a. Prostate covariates (should pass cleanly):\n")
check_covariate_quality(
  prostate_data,
  c("normal_act", "age_cat", "cv_hist", "hemo_bin")
)
cat("  [no warning expected]\n\n")

cat("1b. High-cardinality character column:\n")
df_hc <- prostate_data
df_hc$dummy_chr <- as.character(seq_len(nrow(df_hc)))
suppressWarnings(check_covariate_quality(df_hc, c("dummy_chr")))
tryCatch(
  check_covariate_quality(df_hc, c("dummy_chr")),
  warning = function(w) cat("  caught:", conditionMessage(w), "\n")
)

cat("\n1c. Continuous numeric with many unique values (should NOT warn):\n")
df_cont <- prostate_data
set.seed(1)
df_cont$age_num <- rnorm(nrow(df_cont), 70, 10)
tryCatch(
  {
    check_covariate_quality(df_cont, c("age_num"))
    cat("  [no warning - correct]\n")
  },
  warning = function(w) cat("  UNEXPECTED warning:", conditionMessage(w), "\n")
)

cat("\n1d. Unsupported type (Date):\n")
df_bad <- prostate_data
df_bad$bad_date <- as.Date("2020-01-01") + seq_len(nrow(df_bad))
tryCatch(
  check_covariate_quality(df_bad, c("bad_date")),
  warning = function(w) cat("  caught:", conditionMessage(w), "\n")
)

cat("\n1e. Constant covariate:\n")
df_const <- prostate_data
df_const$dummy_const <- 1L
tryCatch(
  check_covariate_quality(df_const, c("dummy_const")),
  warning = function(w) cat("  caught:", conditionMessage(w), "\n")
)


# =============================================================================
# 2. validate_subject_level (new n_intervals/cut_points args, NA collection)
# =============================================================================
cat("\n\n============================================================\n")
cat(" 2. validate_subject_level\n")
cat("============================================================\n\n")

happy_args <- list(
  data = prostate_data,
  id = "id", time = "event_time", event = "event_type",
  treatment = "A",
  covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
  event_y = 1, event_d = 2, event_c = 0,
  time_varying = NULL
)

cat("2a. Happy path (returns TRUE invisibly):\n")
cat("  result =", do.call(validate_subject_level, happy_args), "\n\n")

cat("2b. NULL id column:\n")
args <- happy_args; args$id <- NULL
tryCatch(
  do.call(validate_subject_level, args),
  error = function(e) cat("  caught:", conditionMessage(e), "\n")
)

cat("\n2c. Missing column:\n")
args <- happy_args; args$time <- "nonexistent"
tryCatch(
  do.call(validate_subject_level, args),
  error = function(e) cat("  caught:", conditionMessage(e), "\n")
)

cat("\n2d. Multiple NAs in critical columns (should report all at once):\n")
df_nas <- prostate_data
df_nas$A[1] <- NA
df_nas$event_time[2] <- NA
args <- happy_args; args$data <- df_nas
tryCatch(
  do.call(validate_subject_level, args),
  error = function(e) cat("  caught:\n", conditionMessage(e), "\n")
)

cat("\n2e. NA-as-censoring pattern (targeted fix message):\n")
df_naev <- prostate_data
df_naev$event_type[df_naev$event_type == 0] <- NA
args <- happy_args; args$data <- df_naev
tryCatch(
  do.call(validate_subject_level, args),
  error = function(e) cat("  caught:\n", conditionMessage(e), "\n")
)

cat("\n2f. Duplicate event labels:\n")
args <- happy_args; args$event_d <- 1  # same as event_y
tryCatch(
  do.call(validate_subject_level, args),
  error = function(e) cat("  caught:", conditionMessage(e), "\n")
)

cat("\n2g. n_intervals + cut_points both given:\n")
args <- happy_args; args$n_intervals <- 12; args$cut_points <- c(10, 20)
tryCatch(
  do.call(validate_subject_level, args),
  error = function(e) cat("  caught:", conditionMessage(e), "\n")
)

cat("\n2h. n_intervals not a positive integer:\n")
args <- happy_args; args$n_intervals <- -1
tryCatch(
  do.call(validate_subject_level, args),
  error = function(e) cat("  caught:", conditionMessage(e), "\n")
)

cat("\n2i. cut_points outside range:\n")
args <- happy_args; args$cut_points <- c(10, 20, 10000)
tryCatch(
  do.call(validate_subject_level, args),
  error = function(e) cat("  caught:", conditionMessage(e), "\n")
)


# =============================================================================
# 3. to_person_time with n_intervals / cut_points
# =============================================================================
cat("\n\n============================================================\n")
cat(" 3. to_person_time (new n_intervals / cut_points)\n")
cat("============================================================\n\n")

cat("3a. Default (n_intervals = 12):\n")
pt_default <- to_person_time(
  prostate_data,
  id = "id", time = "event_time", event = "event_type",
  treatment = "A",
  covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
  event_y = 1, event_d = 2, event_c = 0
)
cat("  class:", paste(class(pt_default), collapse = ", "), "\n")
cat("  rows:", nrow(pt_default), "\n")
cat("  n cut_times:", length(attr(pt_default, "cut_times")), "\n")
cat("  cut_times:", round(attr(pt_default, "cut_times"), 2), "\n\n")

cat("3b. n_intervals = 6:\n")
pt_6 <- to_person_time(
  prostate_data,
  id = "id", time = "event_time", event = "event_type",
  treatment = "A",
  covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
  event_y = 1, event_d = 2, event_c = 0,
  n_intervals = 6
)
cat("  n cut_times:", length(attr(pt_6, "cut_times")), "\n")
cat("  cut_times:", round(attr(pt_6, "cut_times"), 2), "\n\n")

cat("3c. Explicit cut_points = c(12, 24, 36, 48, 60, 72):\n")
pt_cp <- to_person_time(
  prostate_data,
  id = "id", time = "event_time", event = "event_type",
  treatment = "A",
  covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
  event_y = 1, event_d = 2, event_c = 0,
  cut_points = c(12, 24, 36, 48, 60, 72)
)
cat("  cut_times:", attr(pt_cp, "cut_times"), "\n\n")

cat("3d. Event flag integrity (subject 3 has event_type=2, event_time=41):\n")
sub3 <- pt_default[pt_default$id == 3, c("k", "event_type", "y_flag", "d_flag", "c_flag")]
print(sub3)


# =============================================================================
# 4. separable_effects end-to-end (method = "all")
# =============================================================================
cat("\n\n============================================================\n")
cat(" 4. separable_effects (method = 'all')\n")
cat("============================================================\n\n")

fit_all <- separable_effects(pt_default)

cat("4a. class:", paste(class(fit_all), collapse = ", "), "\n")
cat("4b. method:", fit_all$method, "\n")
cat("4c. n:", fit_all$n, "\n")
cat("4d. cumulative_incidence is a list:", is.list(fit_all$cumulative_incidence), "\n")
cat("4e. methods present in $cumulative_incidence:\n")
cat("   names:", names(fit_all$cumulative_incidence), "\n")

cat("\n4f. g-formula cumulative incidence at final time:\n")
last_gf <- fit_all$cumulative_incidence$gformula
last_gf <- last_gf[nrow(last_gf), ]
cat("   arm(1,1) =", round(last_gf$arm_11, 4), "\n")
cat("   arm(0,0) =", round(last_gf$arm_00, 4), "\n")
cat("   arm(1,0) =", round(last_gf$arm_10, 4), "\n")

cat("\n4g. IPW Rep 1 cumulative incidence at final time:\n")
last_ipw1 <- fit_all$cumulative_incidence$ipw1
last_ipw1 <- last_ipw1[nrow(last_ipw1), ]
cat("   arm(1,1) =", round(last_ipw1$arm_11, 4), "\n")
cat("   arm(0,0) =", round(last_ipw1$arm_00, 4), "\n")
cat("   arm(1,0) =", round(last_ipw1$arm_10, 4), "\n")

cat("\n4h. IPW Rep 2 cumulative incidence at final time:\n")
last_ipw2 <- fit_all$cumulative_incidence$ipw2
last_ipw2 <- last_ipw2[nrow(last_ipw2), ]
cat("   arm(1,1) =", round(last_ipw2$arm_11, 4), "\n")
cat("   arm(0,0) =", round(last_ipw2$arm_00, 4), "\n")
cat("   arm(1,0) =", round(last_ipw2$arm_10, 4), "\n")

cat("\n4i. weights slot present:", !is.null(fit_all$weights), "\n")
if (!is.null(fit_all$weights)) {
  cat("   pt_data_weighted cols added:\n")
  base_cols <- names(pt_default)
  new_cols <- setdiff(names(fit_all$weights$pt_data_weighted), base_cols)
  cat("  ", paste(new_cols, collapse = ", "), "\n")
}

cat("\n4j. warnings slot:\n")
if (length(fit_all$warnings) > 0) {
  for (w in fit_all$warnings) cat("   -", w, "\n")
} else {
  cat("   (none)\n")
}


# =============================================================================
# 5. method = "gformula" (no IPW, no weights)
# =============================================================================
cat("\n\n============================================================\n")
cat(" 5. separable_effects (method = 'gformula')\n")
cat("============================================================\n\n")

fit_gf <- separable_effects(pt_default, method = "gformula")
cat("5a. methods present in $cumulative_incidence: ",
    paste(names(fit_gf$cumulative_incidence), collapse = ", "), "\n")
cat("5b. $weights is NULL:", is.null(fit_gf$weights), "\n")


# =============================================================================
# 6. method = "ipw2" only
# =============================================================================
cat("\n\n============================================================\n")
cat(" 6. separable_effects (method = 'ipw2')\n")
cat("============================================================\n\n")

fit_ipw2 <- separable_effects(pt_default, method = "ipw2")
cat("6a. methods present:", paste(names(fit_ipw2$cumulative_incidence), collapse = ", "), "\n")

last_ipw2_only <- fit_ipw2$cumulative_incidence$ipw2
last_ipw2_only <- last_ipw2_only[nrow(last_ipw2_only), ]
cat("6b. Rep 2 cumulative incidence (final time):\n")
cat("   arm(1,1) =", round(last_ipw2_only$arm_11, 4), "\n")
cat("   arm(0,0) =", round(last_ipw2_only$arm_00, 4), "\n")
cat("   arm(1,0) =", round(last_ipw2_only$arm_10, 4), "\n")


# =============================================================================
# 7. censoring_weights = FALSE (skip censoring model)
# =============================================================================
cat("\n\n============================================================\n")
cat(" 7. separable_effects (censoring_weights = FALSE)\n")
cat("============================================================\n\n")

fit_nocw <- separable_effects(pt_default, method = "all", censoring_weights = FALSE)
cat("7a. model_c is NULL:", is.null(fit_nocw$models$model_c), "\n")

if (!is.null(fit_nocw$weights)) {
  w_df <- fit_nocw$weights$pt_data_weighted
  cat("7b. 'w_cens' column present in weighted data:",
      "w_cens" %in% names(w_df), "\n")
  cat("   (should be FALSE - no censoring weights computed)\n")
}

cat("\n7c. IPW Rep 1 final values with censoring_weights = FALSE:\n")
if (!is.null(fit_nocw$cumulative_incidence$ipw1)) {
  last_ipw1_nocw <- fit_nocw$cumulative_incidence$ipw1
  last_ipw1_nocw <- last_ipw1_nocw[nrow(last_ipw1_nocw), ]
  cat("   arm(1,1) =", round(last_ipw1_nocw$arm_11, 4), "\n")
  cat("   arm(0,0) =", round(last_ipw1_nocw$arm_00, 4), "\n")
  cat("   arm(1,0) =", round(last_ipw1_nocw$arm_10, 4), "\n")
}

cat("\n=== All sections complete ===\n")
