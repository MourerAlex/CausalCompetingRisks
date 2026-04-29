#' Validate Basic Input Data Frame Shape
#'
#' Checks that an input object is a non-empty data.frame. Shared by
#' [separable_effects()] and [to_person_time()] as their first line of defense.
#'
#' @param x The object to check.
#' @param name Character. Name of the argument (for error messages).
#'
#' @return Invisibly TRUE if all checks pass; otherwise throws a helpful
#'   error.
#' @family internal
#' @keywords internal
validate_input_shape <- function(x, name) {
  if (is.null(x)) {
    stop(name, " is NULL.", call. = FALSE)
  }
  if (!is.data.frame(x)) {
    stop(
      name, " must be a data.frame. Got class: ",
      paste(class(x), collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(x) == 0) {
    stop(name, " has 0 rows.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Validate Subject-Level Input
#'
#' Checks that subject-level input (one row per subject) and its arguments
#' passed to [to_person_time()] are valid. Raises informative errors for
#' invalid inputs and warnings for potential issues.
#'
#' @param data A data.frame. One row per subject.
#' @param id,time,event,treatment Character column names.
#' @param covariates Character vector of covariate column names.
#' @param event_y,event_d,event_c Values in `event` corresponding to primary
#'   event, competing event, and censoring.
#' @param time_varying Character vector or NULL.
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @details
#' ## Hard errors
#' - NULL column name arguments
#' - Column doesn't exist in data
#' - NAs in critical columns (id, time, event, treatment)
#' - Event column with 2 non-NA values + NAs (likely NA = censored)
#' - Treatment not exactly 2 unique values
#' - Event column not exactly 3 unique values
#' - event_y, event_d, event_c NULL or not matching event column values
#' - Negative time values
#' - time_varying not NULL (not implemented in v1)
#' - Not subject-level (more rows than unique IDs)
#'
#' ## Warnings
#' - Covariates with NAs
#' - Covariate with > 20 unique values
#' - Constant covariates (< 2 unique values)
#' - Event category with < 1% of observations
#' - Subjects with event_time = 0
#'
#' @family internal
#' @keywords internal
validate_subject_level <- function(data,
                                   id,
                                   time,
                                   event,
                                   treatment,
                                   covariates,
                                   event_y,
                                   event_d,
                                   event_c,
                                   n_intervals = NULL,
                                   cut_points = NULL,
                                   time_varying) {

  # --- time_varying not implemented ---
  if (!is.null(time_varying)) {
    stop(
      "Time-varying covariates are not implemented in v1. ",
      "Set time_varying = NULL.",
      call. = FALSE
    )
  }

  # --- Column name arguments must not be NULL ---
  if (is.null(id) || is.null(time) || is.null(event) || is.null(treatment)) {
    stop(
      "id, time, event, and treatment must be column names, not NULL.",
      call. = FALSE
    )
  }

  # --- Event labels must be specified ---
  if (is.null(event_y) || is.null(event_d) || is.null(event_c)) {
    stop(
      "event_y, event_d, and event_c must be specified. ",
      "These are the values in '", event, "' that correspond to the primary ",
      "event, competing event, and censoring.",
      call. = FALSE
    )
  }

  # --- Check columns exist in data ---
  required_cols <- c(id, time, event, treatment, covariates)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) not found in data: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # --- Must be subject-level ---
  n_ids <- length(unique(data[[id]]))
  if (n_ids < nrow(data)) {
    stop(
      "Data has ", nrow(data), " rows but only ", n_ids, " unique IDs. ",
      "validate_subject_level() expects one row per subject. ",
      "If data is already in person-time format, pass it directly to separable_effects().",
      call. = FALSE
    )
  }

  # --- n_intervals / cut_points arg validation ---
  # Note: if both are supplied, cut_points takes precedence (with a warning
  # emitted by to_person_time()). Here we just validate each individually.
  if (!is.null(n_intervals)) {
    if (!is.numeric(n_intervals) || length(n_intervals) != 1 ||
        n_intervals < 1 || n_intervals != round(n_intervals)) {
      stop(
        "n_intervals must be a positive integer.",
        call. = FALSE
      )
    }
  }
  if (!is.null(cut_points)) {
    max_time <- max(data[[time]], na.rm = TRUE)
    if (!is.numeric(cut_points) || length(cut_points) < 1 ||
        any(cut_points <= 0) || any(cut_points > max_time)) {
      stop(
        "cut_points must be positive numeric values within (0, ",
        max_time, "].",
        call. = FALSE
      )
    }
  }

  # --- Critical columns must have no NAs (collect all) ---
  na_issues <- character()
  for (col_name in c(id, time, treatment)) {
    n_na <- sum(is.na(data[[col_name]]))
    if (n_na > 0) {
      na_issues <- c(na_issues,
        paste0("  - '", col_name, "': ", n_na, " NA value(s)")
      )
    }
  }

  # --- Event column: specific NA-as-censoring pattern (targeted fix) ---
  n_na_event <- sum(is.na(data[[event]]))
  event_vals_nona <- unique(data[[event]][!is.na(data[[event]])])

  if (n_na_event > 0 && length(event_vals_nona) == 2) {
    stop(
      "Event column '", event, "' has ", n_na_event, " NAs and 2 non-NA values. ",
      "If NA represents censoring, recode it before calling separable_effects():\n",
      "  data$", event, "[is.na(data$", event, ")] <- 0",
      call. = FALSE
    )
  }
  if (n_na_event > 0) {
    na_issues <- c(na_issues,
      paste0("  - '", event, "': ", n_na_event, " NA value(s)")
    )
  }

  if (length(na_issues) > 0) {
    stop(
      "Critical columns have missing values:\n",
      paste(na_issues, collapse = "\n"),
      call. = FALSE
    )
  }

  # --- Treatment must have exactly 2 unique values ---
  trt_vals <- unique(data[[treatment]])
  if (length(trt_vals) != 2) {
    stop(
      "Treatment column '", treatment, "' must have exactly 2 unique values. ",
      "Found: ", paste(trt_vals, collapse = ", "),
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

  # --- event_y/d/c must be distinct ---
  supplied_labels <- c(event_y, event_d, event_c)
  if (length(unique(supplied_labels)) != 3) {
    stop(
      "event_y, event_d, and event_c must be three distinct values. ",
      "Got: event_y=", event_y, ", event_d=", event_d, ", event_c=", event_c,
      call. = FALSE
    )
  }

  # --- event_y/d/c must match actual values in the event column ---
  if (!all(supplied_labels %in% event_vals)) {
    missing_labels <- setdiff(supplied_labels, event_vals)
    stop(
      "event_y/event_d/event_c value(s) not found in '", event, "' column: ",
      paste(missing_labels, collapse = ", "),
      call. = FALSE
    )
  }

  # --- No negative times ---
  if (any(data[[time]] < 0)) {
    stop(
      "Time column '", time, "' contains negative values.",
      call. = FALSE
    )
  }

  # --- Warnings ---
  check_covariate_quality(data, covariates)

  # Rare event categories (< 1% of observations)
  event_counts <- table(data[[event]])
  n_total <- nrow(data)
  for (i in seq_along(event_counts)) {
    pct <- event_counts[i] / n_total
    if (pct < 0.01) {
      warning(
        "Event value '", names(event_counts)[i], "' represents only ",
        round(pct * 100, 1), "% of observations (n=", event_counts[i], ").",
        call. = FALSE
      )
    }
  }

  # Subjects with event_time = 0
  n_time_zero <- sum(data[[time]] == 0)
  if (n_time_zero > 0) {
    warning(
      n_time_zero, " subject(s) have ", time, " = 0. ",
      "They will be shifted to interval k = 0 (treated as events in the first month).",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate Person-Time Input
#'
#' Checks that user-supplied person-time data (one row per subject-interval)
#' has the structure required by [separable_effects()]. Only called when the input
#' does NOT already come from [to_person_time()] (i.e., no `"causal_cr_pt"`
#' class).
#'
#' @param pt_data A data.frame in person-time format.
#' @param id,treatment Character column names.
#' @param covariates Character vector.
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @details
#' ## Hard errors
#' - NULL column names
#' - Required columns missing (id, treatment, k, y_flag, d_flag, c_flag)
#' - NAs in id or treatment
#' - Flag columns contain values other than 0, 1, or NA
#' - Duplicate (id, k) pairs
#' - Treatment not binary
#'
#' ## Warnings
#' - Covariate quality checks (NAs, cardinality, constant)
#'
#' @family internal
#' @keywords internal
validate_person_time <- function(pt_data,
                                 id,
                                 treatment,
                                 covariates) {

  # --- Column name arguments must not be NULL ---
  if (is.null(id) || is.null(treatment)) {
    stop(
      "id and treatment must be column names, not NULL.",
      call. = FALSE
    )
  }

  # --- Required person-time columns ---
  required <- c(id, treatment, "k", "y_flag", "d_flag", "c_flag", covariates)
  missing_cols <- setdiff(required, names(pt_data))
  if (length(missing_cols) > 0) {
    stop(
      "Person-time data is missing required column(s): ",
      paste(missing_cols, collapse = ", "), ". ",
      "Use to_person_time() to prepare your data.",
      call. = FALSE
    )
  }

  # --- Critical columns must have no NAs ---
  for (col_name in c(id, treatment)) {
    n_na <- sum(is.na(pt_data[[col_name]]))
    if (n_na > 0) {
      stop(
        "Column '", col_name, "' contains ", n_na, " NA value(s). ",
        "Critical columns must have no missing values.",
        call. = FALSE
      )
    }
  }

  # --- Flag columns must contain only {0, 1, NA} ---
  for (flag_col in c("y_flag", "d_flag", "c_flag")) {
    vals <- unique(pt_data[[flag_col]])
    vals_nona <- vals[!is.na(vals)]
    if (!all(vals_nona %in% c(0L, 1L, 0, 1))) {
      stop(
        "Column '", flag_col, "' must contain only 0, 1, or NA. ",
        "Found: ", paste(vals, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # --- Treatment binary ---
  trt_vals <- unique(pt_data[[treatment]])
  if (length(trt_vals) != 2) {
    stop(
      "Treatment column '", treatment, "' must have exactly 2 unique values. ",
      "Found: ", paste(trt_vals, collapse = ", "),
      call. = FALSE
    )
  }

  # --- Duplicate (id, k) check ---
  dupes <- duplicated(pt_data[, c(id, "k")])
  if (any(dupes)) {
    stop(
      "Duplicate (", id, ", k) pairs detected in person-time data. ",
      "Each subject must have at most one row per interval.",
      call. = FALSE
    )
  }

  # --- Left-truncation check: every subject must have a k = 0 row ---
  # Left-truncated data (subjects first observed at k > 0) is not supported
  # in v1 and will not be supported in future versions. This is a structural
  # assumption of the discrete-time pooled logistic framework.
  min_k_per_subject <- tapply(pt_data$k, pt_data[[id]], min)
  if (any(min_k_per_subject > 0)) {
    n_affected <- sum(min_k_per_subject > 0)
    stop(
      "Left-truncated data is not supported. ",
      n_affected, " subject(s) have no row at k = 0. ",
      "Every subject must be observed from the first interval onward.",
      call. = FALSE
    )
  }

  # --- Warnings ---
  check_covariate_quality(pt_data, covariates)

  invisible(TRUE)
}


#' Shared Covariate Quality Checks
#'
#' Emits warnings for covariate-level issues that are identical for both
#' subject-level and person-time input.
#'
#' @param data A data.frame.
#' @param covariates Character vector of covariate column names.
#' @return Invisibly NULL.
#' @family internal
#' @keywords internal
check_covariate_quality <- function(data, covariates) {

  issues <- character()

  for (cov in covariates) {
    col <- data[[cov]]
    col_class <- class(col)[1]

    # --- Type check ---
    # glm natively handles numeric, integer, factor, character, logical.
    # Anything else (Date, list columns, etc.) is unsupported.
    if (!col_class %in% c("numeric", "integer", "factor",
                          "character", "logical")) {
      issues <- c(issues,
        paste0("  - '", cov, "': unsupported type '", col_class,
               "' (expected numeric, integer, factor, character, or logical)")
      )
      next  # skip remaining checks for invalid types
    }

    # --- NAs ---
    if (any(is.na(col))) {
      issues <- c(issues,
        paste0("  - '", cov, "': contains NAs (will be dropped by glm)")
      )
    }

    # --- Degenerate (constant) ---
    n_unique_nona <- length(unique(col[!is.na(col)]))
    if (n_unique_nona < 2) {
      issues <- c(issues,
        paste0("  - '", cov, "': fewer than 2 unique values")
      )
    }

    # --- High cardinality (only for categorical) ---
    # For numeric/integer columns, many unique values is normal (continuous),
    # so we skip this check entirely.
    if (col_class %in% c("factor", "character")) {
      n_levels <- if (col_class == "factor") {
        length(levels(col))
      } else {
        length(unique(col))
      }
      if (n_levels > 20) {
        issues <- c(issues,
          paste0("  - '", cov, "': ", n_levels, " levels ",
                 "(high cardinality - consider grouping)")
        )
      }
    }
  }

  if (length(issues) > 0) {
    warning(
      "Covariate quality issues:\n",
      paste(issues, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(NULL)
}
