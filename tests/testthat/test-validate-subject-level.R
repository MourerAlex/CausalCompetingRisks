# Tests for validate_subject_level (validate.R)

# ---- Fixture: minimal valid subject-level data ----
make_subj <- function(n = 20) {
  set.seed(42)
  data.frame(
    id         = seq_len(n),
    event_time = sample(1:50, n, replace = TRUE),
    event_type = sample(c(0, 1, 2), n, replace = TRUE),
    A          = sample(c(0L, 1L), n, replace = TRUE),
    cov1       = sample(c(0L, 1L), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

valid_args <- function(df = make_subj()) {
  list(
    data = df,
    id = "id", time = "event_time", event = "event_type",
    treatment = "A",
    covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0,
    time_varying = NULL
  )
}


# ---- Happy path ----

test_that("happy path returns TRUE invisibly", {
  expect_true(do.call(validate_subject_level, valid_args()))
})


# ---- NULL args ----

test_that("NULL id errors", {
  args <- valid_args()
  args["id"] <- list(NULL)
  expect_error(
    do.call(validate_subject_level, args),
    "must be column names, not NULL"
  )
})

test_that("NULL event_y errors with distinct message", {
  args <- valid_args()
  args["event_y"] <- list(NULL)
  expect_error(
    do.call(validate_subject_level, args),
    "event_y, event_d, and event_c must be specified"
  )
})


# ---- Missing columns ----

test_that("missing column errors and names the column", {
  args <- valid_args()
  args$time <- "nonexistent"
  expect_error(
    do.call(validate_subject_level, args),
    "Column\\(s\\) not found in data: nonexistent"
  )
})


# ---- Subject-level format check ----

test_that("person-time data rejected with helpful message", {
  # build fake person-time: id repeated
  df <- make_subj(5)
  df <- rbind(df, df)  # 10 rows, 5 unique IDs
  args <- valid_args(df)
  expect_error(
    do.call(validate_subject_level, args),
    "expects one row per subject"
  )
})


# ---- NAs in critical columns ----

test_that("single NA in critical column errors", {
  df <- make_subj()
  df$A[1] <- NA
  args <- valid_args(df)
  expect_error(
    do.call(validate_subject_level, args),
    "'A': 1 NA value"
  )
})

test_that("multiple NAs across critical columns collected into one error", {
  df <- make_subj()
  df$A[1] <- NA
  df$event_time[2] <- NA
  args <- valid_args(df)

  # Both should appear in the message
  err <- tryCatch(do.call(validate_subject_level, args),
                  error = function(e) conditionMessage(e))
  expect_match(err, "event_time")
  expect_match(err, "'A'")
})


# ---- NA-as-censoring targeted fix ----

test_that("event with NAs + 2 non-NA values gives targeted recode message", {
  df <- make_subj()
  df$event_type[df$event_type == 0] <- NA
  args <- valid_args(df)
  err <- tryCatch(do.call(validate_subject_level, args),
                  error = function(e) conditionMessage(e))
  expect_match(err, "recode it before calling separable_effects")
  expect_match(err, "is\\.na\\(data\\$event_type\\)")
})


# ---- Treatment checks ----

test_that("non-binary treatment errors", {
  df <- make_subj()
  df$A[1] <- 2  # 3 unique values
  args <- valid_args(df)
  expect_error(
    do.call(validate_subject_level, args),
    "must have exactly 2 unique values"
  )
})


# ---- Event value checks ----

test_that("event with != 3 unique values errors", {
  df <- make_subj()
  df$event_type <- ifelse(df$event_type == 2, 1, df$event_type)
  args <- valid_args(df)
  expect_error(
    do.call(validate_subject_level, args),
    "must have exactly 3 unique values"
  )
})

test_that("duplicate event labels error", {
  args <- valid_args()
  args$event_d <- 1  # same as event_y
  expect_error(
    do.call(validate_subject_level, args),
    "must be three distinct values"
  )
})

test_that("event labels not matching event column error", {
  args <- valid_args()
  args$event_c <- 99  # not in event column
  expect_error(
    do.call(validate_subject_level, args),
    "not found in 'event_type' column"
  )
})


# ---- Negative times ----

test_that("negative time errors", {
  df <- make_subj()
  df$event_time[1] <- -5
  args <- valid_args(df)
  expect_error(
    do.call(validate_subject_level, args),
    "contains negative values"
  )
})


# ---- time_varying ----

test_that("time_varying != NULL errors (v1 limitation)", {
  args <- valid_args()
  args$time_varying <- "bp"
  expect_error(
    do.call(validate_subject_level, args),
    "not implemented in v1"
  )
})


# ---- n_intervals / cut_points ----

test_that("n_intervals + cut_points both given does NOT error at validate", {
  # validate_subject_level no longer rejects both; to_person_time()
  # resolves precedence with a warning.
  args <- valid_args()
  args$n_intervals <- 12
  args$cut_points <- c(10, 20)
  expect_true(do.call(validate_subject_level, args))
})

test_that("n_intervals = -1 errors", {
  args <- valid_args()
  args$n_intervals <- -1
  expect_error(
    do.call(validate_subject_level, args),
    "positive integer"
  )
})

test_that("n_intervals = 1.5 errors (not integer)", {
  args <- valid_args()
  args$n_intervals <- 1.5
  expect_error(
    do.call(validate_subject_level, args),
    "positive integer"
  )
})

test_that("cut_points outside range errors", {
  args <- valid_args()
  args$cut_points <- c(10, 20, 9999)
  expect_error(
    do.call(validate_subject_level, args),
    "cut_points must be positive numeric values"
  )
})

test_that("cut_points with zero or negative errors", {
  args <- valid_args()
  args$cut_points <- c(-1, 10)
  expect_error(
    do.call(validate_subject_level, args),
    "cut_points must be positive numeric values"
  )
})


# ---- Warnings ----

test_that("rare event (< 1%) warns", {
  # Build data where one event type is just 1 of 200 obs (0.5%)
  set.seed(1)
  df <- data.frame(
    id = seq_len(200),
    event_time = sample(1:50, 200, replace = TRUE),
    event_type = c(1, rep(c(0, 2), times = c(100, 99))),
    A = sample(c(0L, 1L), 200, replace = TRUE),
    cov1 = sample(c(0L, 1L), 200, replace = TRUE)
  )
  args <- valid_args(df)
  expect_warning(
    do.call(validate_subject_level, args),
    "represents only"
  )
})

test_that("event_time = 0 warns", {
  df <- make_subj()
  df$event_time[1] <- 0
  args <- valid_args(df)
  expect_warning(
    do.call(validate_subject_level, args),
    "event_time = 0"
  )
})
