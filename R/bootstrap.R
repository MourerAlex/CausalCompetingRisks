#' Bootstrap Confidence Intervals for Separable Effects
#'
#' Performs subject-level bootstrap resampling to construct confidence
#' intervals for cumulative incidence curves and causal contrasts.
#'
#' @param data,id,time,event,treatment,covariates,event_y,event_d,event_c,method,times,eval_times,formulas,extreme_weight_adjust,extreme_weight_threshold
#'   Arguments passed through from [causal_cr()].
#' @param n_boot Integer. Number of bootstrap replicates.
#' @param alpha Numeric. Significance level for CIs.
#'
#' @return A list with:
#'   \describe{
#'     \item{replicates}{3D array (n_boot x n_arms x n_times) of cumulative
#'       incidence draws.}
#'     \item{ci_curves}{Data frame with lower/upper CI bands on cumulative
#'       incidence for each arm.}
#'     \item{ci_contrasts}{Data frame with CIs on risk differences and ratios
#'       at eval_times.}
#'   }
#'
#' @details
#' ## Resampling scheme
#' Subject-level: sample IDs with replacement, pull all person-time rows for
#' each sampled ID. Duplicate IDs receive a unique suffix to avoid model
#' confusion.
#'
#' ## Per replicate
#' 1. Resample IDs with replacement
#' 2. Rebuild person-time data for resampled subjects
#' 3. Re-fit all hazard models
#' 4. Re-run g-formula and/or IPW estimation
#' 5. Store cumulative incidence for all arms
#'
#' ## CI construction
#' Percentile method: take `alpha/2` and `1 - alpha/2` quantiles of the
#' bootstrap distribution at each time point.
#'
#' @family internal
#' @keywords internal
bootstrap_causal_cr <- function(data, id, time, event, treatment, covariates,
                                event_y, event_d, event_c,
                                method, times, eval_times, formulas,
                                extreme_weight_adjust,
                                extreme_weight_threshold,
                                n_boot, alpha) {

  unique_ids <- unique(data[[id]])
  n <- length(unique_ids)

  # Storage
  arm_names <- c("arm_11", "arm_00", "arm_10")

  # Determine number of time points from first estimation
  prep_first <- to_person_time(
    data, id, time, event, treatment, covariates,
    event_y, event_d, event_c, times
  )
  n_times <- length(prep_first$cut_times)
  boot_array <- array(
    NA_real_,
    dim = c(n_boot, length(arm_names), n_times),
    dimnames = list(NULL, arm_names, NULL)
  )

  for (b in seq_len(n_boot)) {
    if (b %% 50 == 0 || b == 1) {
      message("Bootstrap replicate ", b, "/", n_boot)
    }

    # Resample IDs with replacement
    sampled_ids <- sample(unique_ids, size = n, replace = TRUE)

    # Build resampled dataset with unique IDs
    boot_data <- do.call(rbind, lapply(seq_along(sampled_ids), function(i) {
      rows <- data[data[[id]] == sampled_ids[i], ]
      rows[[id]] <- paste0(sampled_ids[i], "_", i)
      rows
    }))

    # Re-estimate (wrapped in tryCatch for robustness)
    result <- tryCatch({
      prep <- to_person_time(
        boot_data, id, time, event, treatment, covariates,
        event_y, event_d, event_c, times
      )

      models <- fit_hazard_models(
        prep$person_time, treatment, covariates, method, formulas
      )

      if (method %in% c("both", "gformula")) {
        ci <- gformula_estimate(
          prep$person_time, models, treatment, prep$cut_times
        )
        ci
      } else {
        ipw_res <- ipw_estimate(
          prep$person_time, models, treatment, method, prep$cut_times,
          extreme_weight_adjust, extreme_weight_threshold
        )
        ipw_res$cumulative_incidence
      }
    }, error = function(e) NULL)

    if (!is.null(result)) {
      for (j in seq_along(arm_names)) {
        boot_array[b, j, ] <- result[[arm_names[j]]]
      }
    }
  }

  # --- Percentile CIs ---
  lo <- alpha / 2
  hi <- 1 - alpha / 2

  ci_curves <- data.frame(k = prep_first$cut_times)
  for (arm in arm_names) {
    j <- which(arm_names == arm)
    ci_curves[[paste0(arm, "_lower")]] <- apply(
      boot_array[, j, , drop = FALSE], 3,
      stats::quantile, probs = lo, na.rm = TRUE
    )
    ci_curves[[paste0(arm, "_upper")]] <- apply(
      boot_array[, j, , drop = FALSE], 3,
      stats::quantile, probs = hi, na.rm = TRUE
    )
  }

  # --- Contrast CIs ---
  boot_contrasts <- array(NA_real_, dim = c(n_boot, 6, n_times))
  dimnames(boot_contrasts) <- list(
    NULL,
    c("total_rd", "total_rr", "sep_direct_rd", "sep_direct_rr",
      "sep_indirect_rd", "sep_indirect_rr"),
    NULL
  )

  for (b in seq_len(n_boot)) {
    a11 <- boot_array[b, "arm_11", ]
    a00 <- boot_array[b, "arm_00", ]
    a10 <- boot_array[b, "arm_10", ]
    boot_contrasts[b, "total_rd", ] <- a11 - a00
    boot_contrasts[b, "total_rr", ] <- a11 / a00
    boot_contrasts[b, "sep_direct_rd", ] <- a10 - a00
    boot_contrasts[b, "sep_direct_rr", ] <- a10 / a00
    boot_contrasts[b, "sep_indirect_rd", ] <- a11 - a10
    boot_contrasts[b, "sep_indirect_rr", ] <- a11 / a10
  }

  ci_contrasts <- data.frame(k = prep_first$cut_times)
  for (cname in dimnames(boot_contrasts)[[2]]) {
    j <- which(dimnames(boot_contrasts)[[2]] == cname)
    ci_contrasts[[paste0(cname, "_lower")]] <- apply(
      boot_contrasts[, j, , drop = FALSE], 3,
      stats::quantile, probs = lo, na.rm = TRUE
    )
    ci_contrasts[[paste0(cname, "_upper")]] <- apply(
      boot_contrasts[, j, , drop = FALSE], 3,
      stats::quantile, probs = hi, na.rm = TRUE
    )
  }

  list(
    replicates = boot_array,
    ci_curves = ci_curves,
    ci_contrasts = ci_contrasts
  )
}
