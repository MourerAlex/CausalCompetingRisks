#' Plot Cumulative Incidence Curves
#'
#' Plots cumulative incidence for each treatment arm with optional bootstrap
#' confidence bands. Uses the Okabe-Ito colorblind-safe palette.
#'
#' @param x A `"causal_cr_risk"` object from [risk()].
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object.
#'
#' @family plot
#' @export
plot.causal_cr_risk <- function(x, ...) {
  ci <- x$cumulative_incidence
  boot_ci <- x$bootstrap

  # Reshape to long format
  plot_data <- data.frame(
    k = rep(ci$k, 3),
    cum_inc = c(ci$arm_11, ci$arm_00, ci$arm_10),
    arm = rep(
      c("(1,1) Treated", "(0,0) Control", "(1,0) Separable"),
      each = nrow(ci)
    )
  )

  # Okabe-Ito palette
  arm_colors <- c(
    "(1,1) Treated" = "#000000",
    "(0,0) Control" = "#0072B2",
    "(1,0) Separable" = "#009E73"
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$k, y = .data$cum_inc, color = .data$arm
  )) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = arm_colors) +
    ggplot2::labs(
      x = "Time",
      y = "Cumulative Incidence",
      color = "Arm"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  # Add CI ribbons if bootstrap available
  if (!is.null(boot_ci)) {
    ribbon_data <- data.frame(
      k = rep(boot_ci$k, 3),
      lower = c(boot_ci$arm_11_lower, boot_ci$arm_00_lower, boot_ci$arm_10_lower),
      upper = c(boot_ci$arm_11_upper, boot_ci$arm_00_upper, boot_ci$arm_10_upper),
      arm = rep(
        c("(1,1) Treated", "(0,0) Control", "(1,0) Separable"),
        each = nrow(boot_ci)
      )
    )
    p <- p + ggplot2::geom_ribbon(
      data = ribbon_data,
      ggplot2::aes(
        x = .data$k, ymin = .data$lower, ymax = .data$upper,
        fill = .data$arm
      ),
      alpha = 0.2, inherit.aes = FALSE
    ) +
      ggplot2::scale_fill_manual(values = arm_colors, guide = "none")
  }

  p
}


#' Plot Effect-Over-Time Curves
#'
#' Plots risk differences (or risk ratios) for total, separable direct, and
#' separable indirect effects over time. Uses Okabe-Ito palette:
#' total = black, sep direct (A_Y) = green, sep indirect (A_D) = vermillion.
#'
#' @param x A `"causal_cr_contrast"` object from [contrast()].
#' @param type Character. `"rd"` for risk differences (default) or `"rr"` for
#'   risk ratios.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object.
#'
#' @family plot
#' @export
plot.causal_cr_contrast <- function(x, type = "rd", ...) {
  ct <- x$contrasts
  boot_ci <- x$bootstrap

  if (type == "rd") {
    plot_data <- data.frame(
      k = rep(ct$k, 3),
      effect = c(ct$total_rd, ct$sep_direct_rd, ct$sep_indirect_rd),
      type = rep(
        c("Total", "Sep. direct (A_Y)", "Sep. indirect (A_D)"),
        each = nrow(ct)
      )
    )
    ylab <- "Risk Difference"
    ref_line <- 0
  } else {
    plot_data <- data.frame(
      k = rep(ct$k, 3),
      effect = c(ct$total_rr, ct$sep_direct_rr, ct$sep_indirect_rr),
      type = rep(
        c("Total", "Sep. direct (A_Y)", "Sep. indirect (A_D)"),
        each = nrow(ct)
      )
    )
    ylab <- "Risk Ratio"
    ref_line <- 1
  }

  # Okabe-Ito palette (agreed colors)
  effect_colors <- c(
    "Total" = "#000000",
    "Sep. direct (A_Y)" = "#009E73",
    "Sep. indirect (A_D)" = "#D55E00"
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$k, y = .data$effect, color = .data$type
  )) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(values = effect_colors) +
    ggplot2::labs(
      x = "Time",
      y = ylab,
      color = "Effect"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  # Add CI ribbons if bootstrap available
  if (!is.null(boot_ci)) {
    if (type == "rd") {
      cols_lower <- c("total_rd_lower", "sep_direct_rd_lower", "sep_indirect_rd_lower")
      cols_upper <- c("total_rd_upper", "sep_direct_rd_upper", "sep_indirect_rd_upper")
    } else {
      cols_lower <- c("total_rr_lower", "sep_direct_rr_lower", "sep_indirect_rr_lower")
      cols_upper <- c("total_rr_upper", "sep_direct_rr_upper", "sep_indirect_rr_upper")
    }

    ribbon_data <- data.frame(
      k = rep(boot_ci$k, 3),
      lower = c(boot_ci[[cols_lower[1]]], boot_ci[[cols_lower[2]]], boot_ci[[cols_lower[3]]]),
      upper = c(boot_ci[[cols_upper[1]]], boot_ci[[cols_upper[2]]], boot_ci[[cols_upper[3]]]),
      type = rep(
        c("Total", "Sep. direct (A_Y)", "Sep. indirect (A_D)"),
        each = nrow(boot_ci)
      )
    )
    p <- p + ggplot2::geom_ribbon(
      data = ribbon_data,
      ggplot2::aes(
        x = .data$k, ymin = .data$lower, ymax = .data$upper,
        fill = .data$type
      ),
      alpha = 0.2, inherit.aes = FALSE
    ) +
      ggplot2::scale_fill_manual(values = effect_colors, guide = "none")
  }

  p
}


#' Plot Weight Diagnostics
#'
#' Plots weight distributions for IPW diagnostics.
#'
#' @param x A `"causal_cr_diagnostic"` object from [diagnostic()].
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object, or NULL if no weight information is available.
#'
#' @family plot
#' @export
plot.causal_cr_diagnostic <- function(x, ...) {
  if (is.null(x$weight_summary)) {
    message("No IPW diagnostics available (g-formula only?).")
    return(invisible(NULL))
  }

  # Placeholder — will expand with histogram/density of weights
  # and W_D/W_Y departure from 1 over time
  message("Diagnostic plots not yet implemented.")
  invisible(NULL)
}
