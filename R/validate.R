#' Validate Input Data and Arguments
#'
#' Checks that the input data and arguments to [causal_cr()] are valid.
#' Raises informative errors for invalid inputs and warnings for potential
#' issues.
#'
#' @param data A data.frame.
#' @param id,time,event,treatment Character column names.
#' @param covariates Character vector of covariate column names.
#' @param event_y,event_d,event_c Character or NULL event labels.
#' @param method Character estimation method.
#' @param time_varying Character vector or NULL.
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @details
#' ## Hard errors
#' - Column doesn't exist in data
#' - Treatment not binary (0/1)
#' - Event column not exactly 3 unique values
#' - Character event without all of event_y, event_d, event_c specified
#' - Unknown method
#' - Duplicate (id, time) pairs
#' - Negative time values
#' - time_varying not NULL (not implemented in v1)
#'
#' ## Warnings
#' - Covariates with NAs
#' - Covariate with > 20 unique values (likely miscoded continuous)
#' - Fewer than 5 events of either type
#'
#' @family internal
#' @keywords internal
validate_input <- function(data,
                           id,
                           time,
                           event,
                           treatment,
                           covariates,
                           event_y,
                           event_d,
                           event_c,
                           method,
                           time_varying) {

  # --- time_varying not implemented ---
  if (!is.null(time_varying)) {
    stop(
      "Time-varying covariates are not implemented in v1. ",
      "Set time_varying = NULL.",
      call. = FALSE
    )
  }

  # --- Check columns exist ---
  required_cols <- c(id, time, event, treatment, covariates)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) not found in data: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # --- Treatment must be binary ---
  trt_vals <- sort(unique(data[[treatment]]))
  if (!identical(trt_vals, c(0L, 1L)) && !identical(trt_vals, c(0, 1))) {
    stop(
      "Treatment column '", treatment, "' must be binary (0/1). ",
      "Found values: ", paste(trt_vals, collapse = ", "),
      call. = FALSE
    )
  }

  # --- Event must have exactly 3 unique values ---
  event_vals <- unique(data[[event]])
  if (length(event_vals) != 3) {
    stop(
      "Event column '", event, "' must have exactly 3 unique values ",
      "(censoring, primary event, competing event). ",
      "Found ", length(event_vals), " unique value(s): ",
      paste(event_vals, collapse = ", "),
      call. = FALSE
    )
  }

  # --- Character events require explicit labels ---
  if (is.character(data[[event]])) {
    if (is.null(event_y) || is.null(event_d) || is.null(event_c)) {
      stop(
        "Event column '", event, "' is character-valued. ",
        "You must specify event_y, event_d, and event_c explicitly.",
        call. = FALSE
      )
    }
  }

  # --- Method must be valid ---
  valid_methods <- c("both", "gformula", "ipw1", "ipw2")
  if (!method %in% valid_methods) {
    stop(
      "Unknown method '", method, "'. ",
      "Must be one of: ", paste(valid_methods, collapse = ", "),
      call. = FALSE
    )
  }

  # --- No negative times ---
  if (any(data[[time]] < 0, na.rm = TRUE)) {
    stop(
      "Time column '", time, "' contains negative values.",
      call. = FALSE
    )
  }

  # --- Check for duplicate (id, time) if person-time format ---
  n_ids <- length(unique(data[[id]]))
  if (n_ids < nrow(data)) {
    dupes <- duplicated(data[, c(id, time)])
    if (any(dupes)) {
      stop(
        "Duplicate (id, time) pairs detected in person-time data. ",
        "Each subject must have at most one row per time point.",
        call. = FALSE
      )
    }
  }

  # --- Warnings ---
  # Covariates with NAs
  for (cov in covariates) {
    if (any(is.na(data[[cov]]))) {
      warning(
        "Covariate '", cov, "' contains NA values. ",
        "These will be dropped by glm().",
        call. = FALSE
      )
    }
  }

  # High-cardinality covariates
  for (cov in covariates) {
    n_unique <- length(unique(data[[cov]]))
    if (n_unique > 20) {
      warning(
        "Covariate '", cov, "' has ", n_unique, " unique values. ",
        "Consider categorizing if this is meant to be categorical.",
        call. = FALSE
      )
    }
  }

  # Few events
  if (is.numeric(data[[event]])) {
    n_y <- sum(data[[event]] == 1, na.rm = TRUE)
    n_d <- sum(data[[event]] == 2, na.rm = TRUE)
    if (n_y < 5) {
      warning("Fewer than 5 primary events (Y) detected.", call. = FALSE)
    }
    if (n_d < 5) {
      warning("Fewer than 5 competing events (D) detected.", call. = FALSE)
    }
  }

  invisible(TRUE)
}
