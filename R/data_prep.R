#' Prepare Subject-Level Data for [separable_effects()]
#'
#' Converts subject-level data (one row per subject) into discrete-time
#' person-time format using [survival::survSplit()]. Derives the event-flag
#' columns (`y_flag`, `d_flag`, `c_flag`) and applies the temporal ordering
#' (censoring first, then competing event, then primary event) within each
#' interval.
#'
#' The returned object is a standard data.frame with two extras:
#' an S3 class `"person_time"` and attributes `"cut_times"` and
#' `"event_labels"`. When passed to [separable_effects()], these attributes are
#' reused directly and person-time validation is skipped.
#'
#' @param data A data.frame in subject-level format (one row per subject).
#' @param id Character. Name of the subject identifier column.
#' @param time Character. Name of the event/censoring time column.
#' @param event Character. Name of the event type column. Must have exactly
#'   3 unique values.
#' @param treatment Character. Name of the binary treatment column.
#' @param covariates Character vector. Names of baseline covariate columns.
#' @param event_y Value in `event` indicating the primary event.
#' @param event_d Value in `event` indicating the competing event.
#' @param event_c Value in `event` indicating censoring.
#' @param n_intervals Integer. Number of equally-spaced intervals to use when
#'   no explicit `cut_points` are given. Defaults to 12. Mutually exclusive
#'   with `cut_points`.
#' @param cut_points Numeric vector of explicit time cut points. When `NULL`
#'   (default), `n_intervals` is used. Mutually exclusive with `n_intervals`.
#' @param time_varying Not implemented in v1; must be NULL.
#'
#' @return A data.frame of class `c("person_time", "data.frame")` in
#'   person-time format. Attributes include `"cut_times"` (the time grid)
#'   and `"event_labels"` (the three values mapped to Y, D, and censoring).
#'
#' @details
#' ## survSplit convention
#' We set `event_indicator = 1L` for ALL subjects (including censored). This
#' is critical: censored subjects need their terminal row flagged so that
#' `c_flag` can be derived correctly. Event labels are then matched against
#' the original `event` column on that terminal row to populate `y_flag`,
#' `d_flag`, and `c_flag`.
#'
#' Event times are shifted by +1 so that month 0 represents enrollment
#' (no events possible) and month 1 is the first interval at risk. Subjects
#' with `event_time = 0` land in interval k = 0 (treated as events in the
#' first month).
#'
#' ## Temporal ordering
#' Within each interval, events happen in order C, D, Y:
#' - If censored at k (c_flag = 1): d_flag set to NA.
#' - If censored OR competing event at k: y_flag set to NA.
#'
#' NAs define the risk set downstream — `glm()` automatically excludes NAs
#' from the Y and D hazard models respectively.
#'
#' @examples
#' \dontrun{
#' data(prostate_data)
#' pt <- to_person_time(
#'   prostate_data,
#'   id = "id", time = "event_time", event = "event_type",
#'   treatment = "A",
#'   covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
#'   event_y = 1, event_d = 2, event_c = 0
#' )
#' head(pt)
#' attr(pt, "cut_times")
#' }
#'
#' @seealso [separable_effects()]
#' @export
to_person_time <- function(data,
                           id = "id",
                           time = "event_time",
                           event = "event_type",
                           treatment = "A",
                           covariates = character(),
                           event_y,
                           event_d,
                           event_c,
                           n_intervals = NULL,
                           cut_points = NULL,
                           time_varying = NULL) {

  # --- Basic input shape checks ---
  validate_input_shape(data, "data")

  # --- Validate subject-level input (also checks n_intervals / cut_points) ---
  validate_subject_level(
    data = data,
    id = id,
    time = time,
    event = event,
    treatment = treatment,
    covariates = covariates,
    event_y = event_y,
    event_d = event_d,
    event_c = event_c,
    n_intervals = n_intervals,
    cut_points = cut_points,
    time_varying = time_varying
  )

  # --- Build time grid ---
  max_time <- max(data[[time]])

  if (is.null(n_intervals) && is.null(cut_points)) {
    n_intervals <- 12L  # default
  }

  # cut_points takes precedence if both are supplied
  if (!is.null(cut_points)) {
    if (!is.null(n_intervals)) {
      warning(
        "Both 'n_intervals' and 'cut_points' supplied. ",
        "'cut_points' takes precedence; 'n_intervals' ignored.",
        call. = FALSE
      )
    }
    cut_times <- sort(unique(cut_points))
  } else {
    cut_times <- seq(0, max_time, length.out = n_intervals + 1)[-1]
  }

  # --- Prepare for survSplit ---
  # event_indicator = 1L for ALL subjects (including censored).
  # This ensures survSplit flags terminal rows for censored subjects so
  # that c_flag can be derived.
  df <- data
  df$event_indicator <- 1L
  df$tstart <- 0

  # Shift event times +1 so month 0 = enrollment, month 1 = first at risk.
  df$tstop <- df[[time]] + 1

  # --- survSplit ---
  pt_data <- survival::survSplit(
    data  = df,
    cut   = cut_times,
    start = "tstart",
    end   = "tstop",
    event = "event_indicator"
  )

  # --- Interval index ---
  pt_data$k <- pt_data$tstart

  # --- Derive event flags from event column on terminal rows ---
  pt_data$y_flag <- ifelse(
    pt_data$event_indicator == 1 & pt_data[[event]] == event_y, 1L, 0L
  )
  pt_data$d_flag <- ifelse(
    pt_data$event_indicator == 1 & pt_data[[event]] == event_d, 1L, 0L
  )
  pt_data$c_flag <- ifelse(
    pt_data$event_indicator == 1 & pt_data[[event]] == event_c, 1L, 0L
  )

  # --- Temporal ordering: C_k, D_k, Y_k ---
  # If censored at k, D is undefined (NA).
  # Then Y is undefined if d_flag is NA (censored) OR d_flag == 1 (D occurred).
  pt_data$d_flag <- ifelse(pt_data$c_flag == 1L, NA_integer_, pt_data$d_flag)
  pt_data$y_flag <- ifelse(
    is.na(pt_data$d_flag) | pt_data$d_flag == 1L,
    NA_integer_,
    pt_data$y_flag
  )

  # --- Create A_y and A_d working copies for cross-arm prediction ---
  pt_data$A_y <- pt_data[[treatment]]
  pt_data$A_d <- pt_data[[treatment]]

  # --- Attach metadata and class ---
  attr(pt_data, "cut_times") <- cut_times
  attr(pt_data, "event_labels") <- list(y = event_y, d = event_d, c = event_c)
  attr(pt_data, "id_col") <- id
  attr(pt_data, "treatment_col") <- treatment
  attr(pt_data, "covariates") <- covariates
  class(pt_data) <- c("person_time", class(pt_data))

  pt_data
}
