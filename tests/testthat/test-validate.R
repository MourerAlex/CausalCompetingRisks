test_that("validate_input catches missing columns", {
  df <- data.frame(id = 1:5, time = 1:5, event_type = c(0, 1, 2, 0, 1), A = c(0, 1, 0, 1, 0))
  expect_error(
    validate_input(df, "id", "time", "event_type", "A", "nonexistent",
                   NULL, NULL, NULL, "both", NULL),
    "not found"
  )
})

test_that("validate_input catches non-binary treatment", {
  df <- data.frame(id = 1:5, time = 1:5, event_type = c(0, 1, 2, 0, 1), A = c(0, 1, 2, 1, 0))
  expect_error(
    validate_input(df, "id", "time", "event_type", "A", character(),
                   NULL, NULL, NULL, "both", NULL),
    "binary"
  )
})

test_that("validate_input catches wrong number of event values", {
  df <- data.frame(id = 1:5, time = 1:5, event_type = c(0, 1, 0, 0, 1), A = c(0, 1, 0, 1, 0))
  expect_error(
    validate_input(df, "id", "time", "event_type", "A", character(),
                   NULL, NULL, NULL, "both", NULL),
    "exactly 3"
  )
})

test_that("validate_input catches time_varying in v1", {
  df <- data.frame(id = 1:5, time = 1:5, event_type = c(0, 1, 2, 0, 1), A = c(0, 1, 0, 1, 0))
  expect_error(
    validate_input(df, "id", "time", "event_type", "A", character(),
                   NULL, NULL, NULL, "both", "bp"),
    "not implemented"
  )
})

test_that("validate_input catches unknown method", {
  df <- data.frame(id = 1:5, time = 1:5, event_type = c(0, 1, 2, 0, 1), A = c(0, 1, 0, 1, 0))
  expect_error(
    validate_input(df, "id", "time", "event_type", "A", character(),
                   NULL, NULL, NULL, "magic", NULL),
    "Unknown method"
  )
})

test_that("validate_input catches negative times", {
  df <- data.frame(id = 1:5, time = c(-1, 1, 2, 3, 4), event_type = c(0, 1, 2, 0, 1), A = c(0, 1, 0, 1, 0))
  expect_error(
    validate_input(df, "id", "time", "event_type", "A", character(),
                   NULL, NULL, NULL, "both", NULL),
    "negative"
  )
})

test_that("validate_input requires event labels for character events", {
  df <- data.frame(
    id = 1:5, time = 1:5,
    event_type = c("cens", "death_y", "death_d", "cens", "death_y"),
    A = c(0, 1, 0, 1, 0)
  )
  expect_error(
    validate_input(df, "id", "time", "event_type", "A", character(),
                   NULL, NULL, NULL, "both", NULL),
    "event_y, event_d, and event_c"
  )
})

test_that("validate_input passes with valid data", {
  df <- data.frame(id = 1:5, time = 1:5, event_type = c(0, 1, 2, 0, 1), A = c(0, 1, 0, 1, 0))
  expect_true(
    validate_input(df, "id", "time", "event_type", "A", character(),
                   NULL, NULL, NULL, "both", NULL)
  )
})
