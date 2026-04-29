#' Bootstrap Confidence Intervals for a separable_effects Fit
#'
#' Performs subject-level bootstrap resampling on a [separable_effects()] fit and
#' constructs percentile confidence intervals for cumulative incidence
#' curves. Returns a separate `"separable_effects_bootstrap"` object that pairs
#' with the fit — pass both to plotting and contrast functions to get
#' confidence bands.
#'
#' All estimation settings (method, treatment column, covariates, formulas,
#' weight adjustment, censoring weights) are pulled from the fit so the
#' bootstrap replicates re-run exactly what the point-estimate call did.
#'
#' @param fit A `"separable_effects"` object from [separable_effects()].
#' @param n_boot Integer. Number of bootstrap replicates (default 500).
#' @param alpha Numeric. Two-sided significance level (default 0.05 for
#'   95% CIs).
#'
#' @return An S3 object of class `"separable_effects_bootstrap"`:
#'   \describe{
#'     \item{replicates}{4D array `[replicate, method, arm, time]` of
#'       cumulative incidence estimates per bootstrap draw. Dim names on
#'       the method and arm dimensions.}
#'     \item{ci_curves}{Named list (one entry per method) of data.frames
#'       with columns `k`, `{arm}_lower`, `{arm}_upper` for all 4 arms.
#'       All methods (g-formula, IPW Rep 1, IPW Rep 2) now emit `arm_01`,
#'       enabling Decomposition B sensitivity throughout.}
#'     \item{n_boot, alpha}{The settings used.}
#'     \item{fit_call}{Copy of `fit$call` for provenance.}
#'   }
#'
#' @details
#' ## Resampling scheme
#' Subject-level resampling with replacement. For each replicate, unique
#' IDs are sampled with replacement; the person-time rows for each sampled
#' ID are pulled (via a pre-split lookup) and stitched back together with
#' unique synthetic IDs so that glm treats duplicate draws as distinct
#' subjects.
#'
#' ## Estimation per replicate
#' Calls [fit_separable_effects()] directly (not [separable_effects()]), bypassing the
#' user-facing wrapper's validation and warning capture. Warnings during
#' replicates are suppressed (they would otherwise spam).
#'
#' ## CI construction
#' Percentile method: `alpha/2` and `1 - alpha/2` quantiles of the bootstrap
#' distribution at each time point, computed per arm per method.
#'
#' ## Progress reporting
#' For n_boot > 50: prints every 10 replicates for the first 50, then
#' prints a time estimate for the remaining replicates (based on the first
#' 50), then every 100 after. For n_boot <= 50: prints every 10.
#'
#' @seealso [separable_effects()], [fit_separable_effects()], [contrast()], [risk()]
#'
#' @examples
#' \dontrun{
#' fit <- separable_effects(pt)
#' boot <- bootstrap(fit, n_boot = 500)
#' plot(risk(fit), ci = boot)
#' contrast(fit, method = "gformula", ci = boot)
#' }
#'
#' @export
bootstrap <- function(fit, n_boot = 500, alpha = 0.05) {

  stopifnot(inherits(fit, "separable_effects"))
  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 1 ||
      n_boot != round(n_boot)) {
    stop("n_boot must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1).", call. = FALSE)
  }

  # --- Pull settings from fit ---
  pt_data                  <- fit$person_time
  id_col                   <- fit$id_col
  treatment_col            <- fit$treatment_col
  covariates_vec           <- fit$covariates
  cut_times                <- fit$times
  active_methods           <- fit$active_methods
  formulas                 <- fit$formulas
  ipcw        <- fit$ipcw
  truncate                 <- fit$truncate

  unique_ids <- unique(pt_data[[id_col]])
  n <- length(unique_ids)

  arm_names_vec <- arm_names()
  n_methods <- length(active_methods)
  n_arms    <- length(arm_names_vec)
  n_times   <- length(cut_times)

  # 4D array: [replicate, method, arm, time]
  boot_array <- array(
    NA_real_,
    dim = c(n_boot, n_methods, n_arms, n_times),
    dimnames = list(NULL, active_methods, arm_names_vec, NULL)
  )

  # Split once for O(1) lookup inside the loop
  pt_by_id <- split(pt_data, pt_data[[id_col]])

  tic <- proc.time()[["elapsed"]]
  estimate_announced <- FALSE

  for (b in seq_len(n_boot)) {

    # Progress messages
    if (b == 1 || (b <= 50 && b %% 10 == 0)) {
      message("Bootstrap replicate ", b, "/", n_boot)
    }
    if (!estimate_announced && b == 50 && n_boot > 50) {
      elapsed <- proc.time()[["elapsed"]] - tic
      per_rep <- elapsed / 50
      remaining <- (n_boot - 50) * per_rep
      # Format helper: integer minutes + integer seconds, with min dropped
      # when zero. "3 min 36 sec" / "45 sec".
      fmt_duration <- function(sec) {
        m <- as.integer(floor(sec / 60))
        s <- as.integer(round(sec - 60 * m))
        if (m > 0) {
          sprintf("%d min %02d sec", m, s)
        } else {
          sprintf("%d sec", s)
        }
      }
      message(sprintf(
        "Bootstrap: 50 replicates done in %s. Estimated remaining: %s (%d more replicates).",
        fmt_duration(elapsed), fmt_duration(remaining), n_boot - 50
      ))
      estimate_announced <- TRUE
    }
    if (b > 50 && b %% 100 == 0) {
      message("Bootstrap replicate ", b, "/", n_boot)
    }

    # Resample
    sampled_ids <- sample(unique_ids, size = n, replace = TRUE)
    boot_data <- do.call(rbind, lapply(seq_along(sampled_ids), function(i) {
      rows <- pt_by_id[[as.character(sampled_ids[i])]]
      rows[[id_col]] <- paste0(sampled_ids[i], "_", i)
      rows
    }))

    # Re-fit and re-estimate (inner worker only; suppress warnings)
    result <- tryCatch(
      suppressWarnings(fit_separable_effects(
        pt_data = boot_data,
        id_col = id_col,
        treatment_col = treatment_col,
        covariates_vec = covariates_vec,
        cut_times = cut_times,
        active_methods = active_methods,
        formulas = formulas,
        ipcw = ipcw,
        truncate = truncate
      )),
      error = function(e) NULL
    )

    if (!is.null(result)) {
      for (m in active_methods) {
        ci_df <- result$cumulative_incidence[[m]]
        if (is.null(ci_df)) next
        for (arm in arm_names_vec) {
          if (arm %in% names(ci_df)) {
            boot_array[b, m, arm, ] <- ci_df[[arm]]
          }
        }
      }
    }
  }

  # --- Percentile CIs per method per arm ---
  lo <- alpha / 2
  hi <- 1 - alpha / 2

  ci_curves <- list()
  for (m in active_methods) {
    df <- data.frame(k = cut_times)
    for (arm in arm_names_vec) {
      # boot_array[, m, arm, ] is [n_boot x n_times]
      slice <- boot_array[, m, arm, , drop = FALSE]
      lower <- apply(slice, 4, stats::quantile, probs = lo, na.rm = TRUE)
      upper <- apply(slice, 4, stats::quantile, probs = hi, na.rm = TRUE)
      df[[paste0(arm, "_lower")]] <- unname(lower)
      df[[paste0(arm, "_upper")]] <- unname(upper)
    }
    ci_curves[[m]] <- df
  }

  structure(
    list(
      replicates = boot_array,
      ci_curves  = ci_curves,
      n_boot     = n_boot,
      alpha      = alpha,
      fit_call   = fit$call
    ),
    class = "separable_effects_bootstrap"
  )
}


#' Print a separable_effects_bootstrap Object
#'
#' @param x A `"separable_effects_bootstrap"` object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.separable_effects_bootstrap <- function(x, ...) {
  cat("Bootstrap Confidence Intervals (separable_effects)\n")
  cat("------------------------------------------\n")
  cat("Replicates: ", x$n_boot, "\n", sep = "")
  cat("Significance level: ", x$alpha, " (", (1 - x$alpha) * 100, "% CIs)\n",
      sep = "")
  cat("Methods: ", paste(names(x$ci_curves), collapse = ", "), "\n", sep = "")
  cat("\nUse `contrast(fit, method = '<name>', ci = <this>)` for contrast CIs,\n")
  cat("or `plot(risk(fit), ci = <this>)` for curves with bands.\n")
  invisible(x)
}
