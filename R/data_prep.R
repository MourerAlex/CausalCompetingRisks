#' Prepare Person-Time Data
#'
#' Converts subject-level data to person-time format using [survival::survSplit()],
#' or validates existing person-time data. Derives event indicators (y_event,
#' d_event, c_event) and applies temporal ordering (C before D before Y).
#'
#' @param data A data.frame. Subject-level or person-time format.
#' @param id,time,event,treatment Character column names.
#' @param covariates Character vector of covariate column names.
#' @param event_y,event_d,event_c Character or NULL event labels.
#' @param times Numeric. Time grid specification (NULL, scalar, or vector).
#'
#' @return A list with:
#'   \describe{
#'     \item{person_time}{Data frame in person-time format with columns
#'       `y_event`, `d_event`, `c_event`, `A_y`, `A_d`, and `k` (interval index).}
#'     \item{cut_times}{Numeric vector of time cut points used.}
#'   }
#'
#' @details
#' ## Data format detection
#' - If `n_distinct(id) == nrow(data)`: subject-level, expanded via survSplit.
#' - If `n_distinct(id) < nrow(data)`: person-time, used as-is.
#'
#' ## survSplit convention
#' Uses `event_indicator = 1L` for ALL subjects including censored. This is
#' critical: censored subjects need their terminal row flagged so that
#' `c_event` can be derived correctly.
#'
#' ## Temporal ordering
#' Within each interval, events are ordered C, D, Y:
#' - Censored at k: `d_event` and `y_event` set to NA
#' - Competing event at k: `y_event` set to NA
#'
#' NAs define risk sets for `glm()` via `na.action`.
#'
#' @family internal
#' @keywords internal
to_person_time <- function(data,
                           id,
                           time,
                           event,
                           treatment,
                           covariates,
                           event_y,
                           event_d,
                           event_c,
                           times) {

  # --- Determine event coding ---
  if (is.null(event_y)) {
    # Integer 0/1/2 convention
    message("Using integer event convention: 0 = censored, 1 = Y, 2 = D")
    cens_val <- 0
    y_val <- 1
    d_val <- 2
  } else {
    cens_val <- event_c
    y_val <- event_y
    d_val <- event_d
  }

  # --- Detect data format ---
  n_ids <- length(unique(data[[id]]))
  is_subject_level <- (n_ids == nrow(data))

  if (is_subject_level) {
    # --- Build time grid ---
    max_time <- max(data[[time]], na.rm = TRUE)

    if (is.null(times)) {
      cut_times <- seq(0, max_time, length.out = 13)[-1]
    } else if (length(times) == 1) {
      cut_times <- seq(0, max_time, length.out = times + 1)[-1]
    } else {
      cut_times <- times
    }
    cut_times <- sort(unique(cut_times))

    # --- Prepare for survSplit ---
    # event_indicator = 1L for ALL subjects (including censored)
    # This ensures survSplit flags terminal rows for censored subjects,
    # which is required to derive c_event.
    df <- data
    df$event_indicator <- 1L
    df$tstart <- 0

    # Shift event times +1 so month 0 = enrollment, month 1 = first at risk
    df$tstop <- df[[time]] + 1

    # --- survSplit ---
    pt_data <- survival::survSplit(
      formula = survival::Surv(tstart, tstop, event_indicator) ~ .,
      data = df,
      cut = cut_times
    )

    # --- Interval index ---
    pt_data$k <- pt_data$tstart

    # --- Derive event indicators from event_type + terminal row ---
    pt_data$y_event <- ifelse(
      pt_data$event_indicator == 1 & pt_data[[event]] == y_val, 1L, 0L
    )
    pt_data$d_event <- ifelse(
      pt_data$event_indicator == 1 & pt_data[[event]] == d_val, 1L, 0L
    )
    pt_data$c_event <- ifelse(
      pt_data$event_indicator == 1 & pt_data[[event]] == cens_val, 1L, 0L
    )

    # --- Temporal ordering: C_k, D_k, Y_k ---
    # Censored at k -> D and Y undefined (NA)
    pt_data$d_event <- ifelse(pt_data$c_event == 1L, NA_integer_, pt_data$d_event)
    pt_data$y_event <- ifelse(pt_data$c_event == 1L, NA_integer_, pt_data$y_event)
    # D occurred at k -> Y undefined (NA)
    pt_data$y_event <- ifelse(
      !is.na(pt_data$d_event) & pt_data$d_event == 1L,
      NA_integer_,
      pt_data$y_event
    )

  } else {
    # Person-time data: validate and use as-is
    pt_data <- data
    cut_times <- sort(unique(pt_data[[time]]))

    # Ensure required columns exist
    for (col in c("y_event", "d_event", "c_event")) {
      if (!col %in% names(pt_data)) {
        stop("Person-time data must contain column '", col, "'.", call. = FALSE)
      }
    }

    if (!"k" %in% names(pt_data)) {
      pt_data$k <- pt_data[[time]]
    }
  }

  # --- Create A_y and A_d copies ---
  pt_data$A_y <- pt_data[[treatment]]
  pt_data$A_d <- pt_data[[treatment]]

  list(
    person_time = pt_data,
    cut_times = cut_times
  )
}
