# Tests for check_covariate_quality (validate.R)

make_df <- function() {
  data.frame(
    num_cont  = rnorm(50),
    int_cov   = sample(c(0L, 1L), 50, replace = TRUE),
    fac_cov   = factor(sample(letters[1:3], 50, replace = TRUE)),
    chr_cov   = sample(letters[1:3], 50, replace = TRUE),
    lgl_cov   = sample(c(TRUE, FALSE), 50, replace = TRUE),
    stringsAsFactors = FALSE
  )
}


test_that("clean covariates produce no warning", {
  df <- make_df()
  expect_silent(
    check_covariate_quality(df, c("num_cont", "int_cov", "fac_cov",
                                   "chr_cov", "lgl_cov"))
  )
})


test_that("covariate with NAs warns", {
  df <- make_df()
  df$int_cov[1:3] <- NA
  expect_warning(
    check_covariate_quality(df, "int_cov"),
    "contains NAs"
  )
})


test_that("constant covariate warns", {
  df <- make_df()
  df$const <- 1L
  expect_warning(
    check_covariate_quality(df, "const"),
    "fewer than 2 unique values"
  )
})


test_that("high-cardinality character column warns", {
  df <- make_df()
  df$hc_chr <- as.character(seq_len(nrow(df)))
  expect_warning(
    check_covariate_quality(df, "hc_chr"),
    "high cardinality"
  )
})


test_that("high-cardinality factor (by levels) warns", {
  df <- make_df()
  df$hc_fac <- factor(as.character(seq_len(nrow(df))))
  expect_warning(
    check_covariate_quality(df, "hc_fac"),
    "high cardinality"
  )
})


test_that("continuous numeric with many unique values does NOT warn", {
  df <- make_df()
  # num_cont has ~50 unique values but is numeric -> no warning expected
  expect_silent(
    check_covariate_quality(df, "num_cont")
  )
})


test_that("unsupported type (Date) produces targeted warning", {
  df <- make_df()
  df$bad_date <- as.Date("2020-01-01") + seq_len(nrow(df))
  expect_warning(
    check_covariate_quality(df, "bad_date"),
    "unsupported type 'Date'"
  )
})


test_that("unsupported type skips further checks for that covariate", {
  df <- make_df()
  # Date column with NAs: without the skip, we'd also warn about NAs;
  # with the skip, only the type warning fires.
  df$bad_date <- as.Date("2020-01-01") + seq_len(nrow(df))
  df$bad_date[1:3] <- NA

  w <- tryCatch(
    check_covariate_quality(df, "bad_date"),
    warning = function(cond) conditionMessage(cond)
  )
  expect_match(w, "unsupported type")
  expect_false(grepl("contains NAs", w))
})


test_that("multiple issues across covariates are grouped into one warning", {
  df <- make_df()
  df$int_cov[1:3] <- NA
  df$const <- 1L

  w <- tryCatch(
    check_covariate_quality(df, c("int_cov", "const")),
    warning = function(cond) conditionMessage(cond)
  )
  expect_match(w, "int_cov.*contains NAs")
  expect_match(w, "const.*fewer than 2 unique values")
})
