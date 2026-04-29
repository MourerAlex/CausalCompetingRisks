# ── Load and preprocess ──

library(Hmisc)
library(dplyr)

getHdata(prostate)
prostate <- as.data.frame(prostate)
prostate$rx     <- as.character(prostate$rx)
prostate$status <- as.character(prostate$status)
prostate$pf     <- as.character(prostate$pf)

prostate_data <- prostate %>%
  filter(rx %in% c("placebo", "5.0 mg estrogen")) %>%
  mutate(
    A = as.integer(rx == "5.0 mg estrogen"),
    event_type = case_when(
      status == "dead - prostatic ca" ~ 1L,
      status == "alive"               ~ 0L,
      TRUE                            ~ 2L),
    event_time = dtime + 1,
    normal_act = as.integer(pf == "normal activity"),
    age_cat    = as.character(cut2(age, c(0, 60, 70, 80, 100))),
    cv_hist    = hx,
    hemo_bin   = as.integer(hg < 12)
  ) %>%
  select(id = patno, A, event_time, event_type,
         normal_act, age_cat, cv_hist, hemo_bin)
# Strip all Hmisc labels
for (col in names(prostate_data)) {
  attr(prostate_data[[col]], "label") <- NULL
  class(prostate_data[[col]]) <- setdiff(class(prostate_data[[col]]), "labelled")
}

save(prostate_data, file = "data/prostate_data.rda", compress = "xz")
str(prostate_data)
