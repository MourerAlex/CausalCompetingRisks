# ============================================================
# test_validate.R
# Interactive test script for validate_subject_level() and
# validate_person_time(). Run line by line in RStudio.
# ============================================================

load("data/prostate_data.rda")
source("R/utils.R")
source("R/validate.R")
source("R/data_prep.R")

# ============================================================
# validate_subject_level()
# ============================================================

# --- Should pass cleanly ---
print(validate_subject_level(
  data = prostate_data,
  id = "id", time = "event_time", event = "event_type",
  treatment = "A",
  covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
  event_y = 1, event_d = 2, event_c = 0,
  time_varying = NULL
))

# --- Hard errors ---

# 1. NULL column name
validate_subject_level(prostate_data, id = NULL, time = "event_time",
  event = "event_type", treatment = "A", covariates = character(),
  event_y = 1, event_d = 2, event_c = 0, time_varying = NULL)

# 2. Missing event labels
validate_subject_level(prostate_data, "id", "event_time", "event_type", "A",
  character(), event_y = NULL, event_d = 2, event_c = 0, time_varying = NULL)

# 3. Missing column
validate_subject_level(prostate_data, "id", "blah", "event_type", "A",
  character(), 1, 2, 0, NULL)

# 4. NA in critical column
df_na <- prostate_data
df_na$A[1] <- NA
validate_subject_level(df_na, "id", "event_time", "event_type", "A",
  character(), 1, 2, 0, NULL)

# 5. NA-as-censoring pattern
df_nacens <- prostate_data
df_nacens$event_type[df_nacens$event_type == 0] <- NA
validate_subject_level(df_nacens, "id", "event_time", "event_type", "A",
  character(), 1, 2, 0, NULL)

# 6. Treatment with 3 values
df_trt <- prostate_data
df_trt$A[1] <- 2
validate_subject_level(df_trt, "id", "event_time", "event_type", "A",
  character(), 1, 2, 0, NULL)

# 7. Event label does not appear in column
validate_subject_level(prostate_data, "id", "event_time", "event_type", "A",
  character(), event_y = 99, event_d = 2, event_c = 0, time_varying = NULL)

# 8. Negative time
df_neg <- prostate_data
df_neg$event_time[1] <- -5
validate_subject_level(df_neg, "id", "event_time", "event_type", "A",
  character(), 1, 2, 0, NULL)

# 9. time_varying set
validate_subject_level(prostate_data, "id", "event_time", "event_type", "A",
  character(), 1, 2, 0, "bp")

# 10. Person-time data passed accidentally
df_pt <- rbind(prostate_data, prostate_data[1, ])
validate_subject_level(df_pt, "id", "event_time", "event_type", "A",
  character(), 1, 2, 0, NULL)

# --- Warnings ---

# 11. Covariate with NAs
df_wna <- prostate_data
df_wna$cv_hist[1:3] <- NA
validate_subject_level(df_wna, "id", "event_time", "event_type", "A",
  c("cv_hist"), 1, 2, 0, NULL)

# 12. Constant covariate
df_const <- prostate_data
df_const$dummy <- 1
validate_subject_level(df_const, "id", "event_time", "event_type", "A",
  c("dummy"), 1, 2, 0, NULL)

# 13. Subjects with event_time = 0
df_zero <- prostate_data
df_zero$event_time[1:5] <- 0
validate_subject_level(df_zero, "id", "event_time", "event_type", "A",
  character(), 1, 2, 0, NULL)


# ============================================================
# validate_person_time()
# ============================================================

pt <- to_person_time(
  prostate_data, "id", "event_time", "event_type", "A",
  c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
  event_y = 1, event_d = 2, event_c = 0
)

# 14. Clean person-time
print(validate_person_time(pt, "id", "A",
  c("normal_act", "age_cat", "cv_hist", "hemo_bin")))

# 15. Missing required column
pt_bad <- pt
pt_bad$y_flag <- NULL
validate_person_time(pt_bad, "id", "A", character())

# 16. Invalid flag value
pt_bad2 <- pt
pt_bad2$y_flag[1] <- 5
validate_person_time(pt_bad2, "id", "A", character())
