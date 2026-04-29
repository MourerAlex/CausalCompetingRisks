# Tests for validate_person_time (validate.R)

# ---- Fixture: minimal valid person-time data ----
make_pt <- function() {
  data.frame(
    id       = rep(1:3, each = 3),
    A        = rep(c(1L, 0L, 1L), each = 3),
    k        = rep(0:2, times = 3),
    y_flag   = c(0L, 0L, 1L,   0L, 1L, NA, 0L, 0L, 0L),
    d_flag   = c(0L, 0L, 0L,   0L, 0L, NA, 0L, 0L, 1L),
    c_flag   = c(0L, 0L, 0L,   0L, 0L, 1L, 0L, 0L, 0L),
    cov1     = rep(c(1L, 0L, 1L), each = 3),
    stringsAsFactors = FALSE
  )
}


# ---- Happy path ----

test_that("happy path returns TRUE invisibly", {
  pt <- make_pt()
  expect_true(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = "cov1")
  )
})


# ---- NULL args ----

test_that("NULL id errors", {
  pt <- make_pt()
  expect_error(
    validate_person_time(pt, id = NULL, treatment = "A",
                         covariates = character()),
    "must be column names, not NULL"
  )
})

test_that("NULL treatment errors", {
  pt <- make_pt()
  expect_error(
    validate_person_time(pt, id = "id", treatment = NULL,
                         covariates = character()),
    "must be column names, not NULL"
  )
})


# ---- Missing required columns ----

test_that("missing k column errors and names it", {
  pt <- make_pt()
  pt$k <- NULL
  expect_error(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = character()),
    "missing required column\\(s\\): k"
  )
})

test_that("missing flag columns all reported", {
  pt <- make_pt()
  pt$y_flag <- NULL
  pt$d_flag <- NULL
  err <- tryCatch(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = character()),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "y_flag")
  expect_match(err, "d_flag")
})

test_that("missing covariate column errors", {
  pt <- make_pt()
  expect_error(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = c("cov1", "nonexistent")),
    "missing required column\\(s\\): nonexistent"
  )
})


# ---- NAs in critical columns ----

test_that("NA in id errors", {
  pt <- make_pt()
  pt$id[1] <- NA
  expect_error(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = character()),
    "Column 'id' contains 1 NA value"
  )
})

test_that("NA in treatment errors", {
  pt <- make_pt()
  pt$A[c(2, 5)] <- NA
  expect_error(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = character()),
    "Column 'A' contains 2 NA value"
  )
})


# ---- Flag value integrity ----

test_that("y_flag with value 2 errors", {
  pt <- make_pt()
  pt$y_flag[1] <- 2L
  expect_error(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = character()),
    "Column 'y_flag' must contain only 0, 1, or NA"
  )
})

test_that("c_flag with character value errors", {
  pt <- make_pt()
  pt$c_flag[1] <- "bad"
  expect_error(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = character()),
    "Column 'c_flag' must contain only 0, 1, or NA"
  )
})


# ---- Treatment checks ----

test_that("non-binary treatment errors", {
  pt <- make_pt()
  pt$A[1] <- 2L
  expect_error(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = character()),
    "must have exactly 2 unique values"
  )
})


# ---- Duplicate (id, k) ----

test_that("duplicate (id, k) pairs error", {
  pt <- make_pt()
  pt <- rbind(pt, pt[1, ])
  expect_error(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = character()),
    "Duplicate \\(id, k\\) pairs"
  )
})


# ---- Left-truncation rejection ----

test_that("subject missing k = 0 row errors", {
  pt <- make_pt()
  # Remove subject 2's k = 0 row (left-truncated)
  pt <- pt[!(pt$id == 2 & pt$k == 0), ]
  expect_error(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = character()),
    "Left-truncated data is not supported"
  )
})

test_that("multiple left-truncated subjects are counted", {
  pt <- make_pt()
  # Remove k = 0 for subjects 1 and 2
  pt <- pt[!(pt$id %in% c(1, 2) & pt$k == 0), ]
  err <- tryCatch(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = character()),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "2 subject\\(s\\) have no row at k = 0")
})


# ---- Covariate quality warning ----

test_that("constant covariate triggers quality warning", {
  pt <- make_pt()
  pt$cov_const <- 1L
  expect_warning(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = c("cov1", "cov_const")),
    "Covariate quality issues"
  )
})

test_that("clean covariates produce no warning", {
  pt <- make_pt()
  expect_silent(
    validate_person_time(pt, id = "id", treatment = "A",
                         covariates = "cov1")
  )
})
