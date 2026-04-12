#' Comprehensive surveillance quality dashboard
#'
#' Produces a multi-panel display combining calibration diagnostics,
#' detection power, estimation quality, and current variant landscape
#' into a single figure suitable for weekly surveillance reports.
#' Designed for programme managers rather than statisticians.
#'
#' @param fit An \code{lfq_fit} object from \code{\link{fit_model}}.
#' @param data The \code{lfq_data} object used for fitting.
#' @param bt Optional \code{lfq_backtest} object for calibration
#'   panel. If \code{NULL}, the calibration panel is omitted.
#' @param target_prevalence Prevalence for detection power
#'   calculation. Default 0.01 (1\%).
#'
#' @return A list of ggplot objects with class
#'   \code{surveillance_dashboard}. A print method renders all
#'   panels.
#'
#' @details
#' The dashboard contains up to four panels:
#' \enumerate{
#'   \item \strong{Current landscape}: Frequency trajectory from
#'     the fitted model.
#'   \item \strong{Growth advantages}: Forest plot of relative
#'     fitness.
#'   \item \strong{Detection power}: Probability of detecting a
#'     1\% variant as a function of sample size.
#'   \item \strong{Calibration}: Reliability diagram (if backtest
#'     data provided).
#' }
#'
#' @seealso \code{\link{surveillance_value}} for EVOI analysis,
#'   \code{\link{alert_threshold}} for sequential detection.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.3, "B" = 0.9),
#'   n_timepoints = 15, seed = 1)
#' fit <- fit_model(sim, engine = "mlr")
#' panels <- surveillance_dashboard(fit, sim)
#' }
#'
#' @export
surveillance_dashboard <- function(fit, data,
                                   bt = NULL,
                                   target_prevalence = 0.01) {

  if (!inherits(fit, "lfq_fit")) {
    cli::cli_abort("{.arg fit} must be an {.cls lfq_fit} object.")
  }

  panels <- list()

  # Panel 1: Current frequency landscape
  panels$landscape <- autoplot(fit, type = "frequency")

  # Panel 2: Growth advantages
  panels$advantages <- tryCatch(
    autoplot(fit, type = "advantage"),
    error = function(e) NULL
  )

  # Panel 3: Detection power curve
  prev_range <- seq(0.005, 0.10, by = 0.005)
  n_total <- fit$nobs / fit$n_timepoints  # avg sequences per period

  power_data <- tibble::tibble(
    prevalence = prev_range,
    detection_prob = stats::pbinom(
      0, size = as.integer(n_total), prob = prev_range,
      lower.tail = FALSE
    )
  )

  panels$detection <- ggplot2::ggplot(power_data,
    ggplot2::aes(x = .data$prevalence, y = .data$detection_prob)) +
    ggplot2::geom_line(colour = "#4daf4a", linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_vline(xintercept = target_prevalence,
                        linetype = "dotted", colour = "#d6604d") +
    ggplot2::scale_x_continuous(labels = function(x) paste0(x * 100, "%")) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Variant prevalence",
      y = "Detection probability",
      title = "Detection power"
    ) +
    ggplot2::theme_minimal(base_size = 10)

  # Panel 4: Calibration (if backtest available)
  if (!is.null(bt) && inherits(bt, "lfq_backtest") && nrow(bt) > 0) {
    cal <- calibrate(bt)
    panels$calibration <- plot(cal, type = "reliability")
  }

  structure(panels, class = "surveillance_dashboard")
}


#' @export
print.surveillance_dashboard <- function(x, ...) {
  cli::cli_h3("Surveillance dashboard")
  cli::cli_text("{length(x)} panels: {paste(names(x), collapse = ', ')}")
  # Display each panel
  for (nm in names(x)) {
    if (!is.null(x[[nm]])) {
      print(x[[nm]])
    }
  }
  invisible(x)
}
