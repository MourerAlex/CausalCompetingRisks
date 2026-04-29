# Tests for to_person_time (data_prep.R)

# ---- Fixture: minimal valid subject-level data ----
make_subj <- function(n = 20, seed = 42) {
  set.seed(seed)
  data.frame(
    id         = seq_len(n),
    event_time = sample(1:50, n, replace = TRUE),
    event_type = sample(c(0L, 1L, 2L), n, replace = TRUE),
    A          = sample(c(0L, 1L), n, replace = TRUE),
    cov1       = sample(c(0L, 1L), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}


test_that("returns object of class causal_cr_pt", {
  pt <- to_person_time(
    make_subj(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0
  )
  expect_s3_class(pt, "causal_cr_pt")
  expect_s3_class(pt, "data.frame")
})


test_that("attaches cut_times attribute", {
  pt <- to_person_time(
    make_subj(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0
  )
  ct <- attr(pt, "cut_times")
  expect_type(ct, "double")
  expect_true(length(ct) > 0)
  expect_true(all(diff(ct) > 0))  # sorted ascending
})


test_that("attaches event_labels attribute", {
  pt <- to_person_time(
    make_subj(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0
  )
  labs <- attr(pt, "event_labels")
  expect_equal(labs$y, 1)
  expect_equal(labs$d, 2)
  expect_equal(labs$c, 0)
})


test_that("attaches id/treatment/covariates attributes", {
  pt <- to_person_time(
    make_subj(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = c("cov1"),
    event_y = 1, event_d = 2, event_c = 0
  )
  expect_equal(attr(pt, "id_col"), "id")
  expect_equal(attr(pt, "treatment_col"), "A")
  expect_equal(attr(pt, "covariates"), "cov1")
})


# ---- Time grid options ----

test_that("default n_intervals = 12", {
  pt <- to_person_time(
    make_subj(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0
  )
  expect_length(attr(pt, "cut_times"), 12L)
})

test_that("n_intervals = 6 produces 6 cut_times", {
  pt <- to_person_time(
    make_subj(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0,
    n_intervals = 6
  )
  expect_length(attr(pt, "cut_times"), 6L)
})

test_that("explicit cut_points are used and sorted", {
  pt <- to_person_time(
    make_subj(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0,
    cut_points = c(20, 5, 10)  # unsorted
  )
  expect_equal(attr(pt, "cut_times"), c(5, 10, 20))
})


# ---- Required columns + working copies ----

test_that("adds required person-time columns", {
  pt <- to_person_time(
    make_subj(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0
  )
  for (col in c("k", "y_flag", "d_flag", "c_flag", "A_y", "A_d")) {
    expect_true(col %in% names(pt),
                info = paste("missing column:", col))
  }
})

test_that("A_y and A_d initially equal treatment", {
  pt <- to_person_time(
    make_subj(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0
  )
  expect_identical(pt$A_y, pt$A)
  expect_identical(pt$A_d, pt$A)
})


# ---- Event flag integrity ----

# Helper: minimal 3-subject fixture (one per event type) to keep treatment
# binary and all three event labels represented.
make_triplet <- function() {
  data.frame(
    id         = c(1L, 2L, 3L),
    event_time = c(10L, 10L, 10L),
    event_type = c(1L, 2L, 0L),   # Y, D, censored
    A          = c(1L, 0L, 1L),
    cov1       = c(1L, 0L, 1L),
    stringsAsFactors = FALSE
  )
}


test_that("Y event subject gets y_flag=1 on terminal row; 0 before", {
  pt <- to_person_time(
    make_triplet(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0,
    n_intervals = 3
  )
  sub <- pt[pt$id == 1, ]  # the Y-event subject
  terminal <- sub[nrow(sub), ]
  expect_equal(terminal$y_flag, 1L)
  expect_equal(terminal$d_flag, 0L)
  expect_equal(terminal$c_flag, 0L)
  non_terminal <- sub[-nrow(sub), ]
  expect_true(all(non_terminal$y_flag == 0L))
  expect_true(all(non_terminal$d_flag == 0L))
  expect_true(all(non_terminal$c_flag == 0L))
})


test_that("D event terminal row has y_flag NA (temporal ordering)", {
  pt <- to_person_time(
    make_triplet(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0,
    n_intervals = 3
  )
  sub <- pt[pt$id == 2, ]  # the D-event subject
  terminal <- sub[nrow(sub), ]
  expect_true(is.na(terminal$y_flag))
  expect_equal(terminal$d_flag, 1L)
  expect_equal(terminal$c_flag, 0L)
})


test_that("censored subject terminal row has y_flag and d_flag NA", {
  pt <- to_person_time(
    make_triplet(),
    id = "id", time = "event_time", event = "event_type",
    treatment = "A", covariates = "cov1",
    event_y = 1, event_d = 2, event_c = 0,
    n_intervals = 3
  )
  sub <- pt[pt$id == 3, ]  # the censored subject
  terminal <- sub[nrow(sub), ]
  expect_true(is.na(terminal$y_flag))
  expect_true(is.na(terminal$d_flag))
  expect_equal(terminal$c_flag, 1L)
})


# ---- Validation propagation ----

test_that("invalid subject-level input bubbles up from validate", {
  df <- make_subj()
  df$A[1] <- 2  # non-binary
  expect_error(
    to_person_time(
      df,
      id = "id", time = "event_time", event = "event_type",
      treatment = "A", covariates = "cov1",
      event_y = 1, event_d = 2, event_c = 0
    ),
    "must have exactly 2 unique values"
  )
})


# ---- event_time = 0 shift ----

test_that("event_time = 0 works (warning + k = 0 terminal row)", {
  # Triplet: subject 1 has event_time = 0; others keep grid non-degenerate
  df <- data.frame(
    id         = c(1L, 2L, 3L),
    event_time = c(0L, 10L, 10L),
    event_type = c(1L, 2L, 0L),
    A          = c(1L, 0L, 1L),
    cov1       = c(1L, 0L, 1L),
    stringsAsFactors = FALSE
  )
  expect_warning(
    pt <- to_person_time(
      df,
      id = "id", time = "event_time", event = "event_type",
      treatment = "A", covariates = "cov1",
      event_y = 1, event_d = 2, event_c = 0,
      n_intervals = 3
    ),
    "event_time = 0"
  )
  sub <- pt[pt$id == 1, ]
  # Subject 1 should have exactly one row at k = 0 (terminal) with y_flag = 1
  expect_equal(nrow(sub), 1L)
  expect_equal(sub$k, 0)
  expect_equal(sub$y_flag, 1L)
})
