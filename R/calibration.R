#' Calibration diagnostics for lineage frequency forecasts
#'
#' Assesses whether prediction intervals from backtesting are
#' well-calibrated by computing PIT (Probability Integral Transform)
#' values, reliability diagrams, and uniformity tests. A perfectly
#' calibrated forecaster produces PIT values that are uniformly
#' distributed on \eqn{[0, 1]}.
#'
#' @param x An \code{lfq_backtest} object from \code{\link{backtest}},
#'   or an \code{lfq_forecast} object from \code{\link{forecast}}
#'   together with observed data.
#' @param observed For \code{lfq_forecast} objects: a numeric vector
#'   of observed frequencies corresponding to the forecast rows.
#'   Ignored for \code{lfq_backtest} objects (which carry their own
#'   observed values).
#' @param n_bins Number of bins for the PIT histogram and reliability
#'   diagram. Default 10.
#'
#' @return A \code{calibration_report} object (S3 class) with
#'   components:
#'   \describe{
#'     \item{pit_values}{Numeric vector of PIT values in \eqn{[0,1]}.}
#'     \item{pit_histogram}{Tibble with \code{bin}, \code{count},
#'       \code{density}, \code{expected} columns.}
#'     \item{reliability}{Tibble with \code{nominal}, \code{observed}
#'       coverage at each level.}
#'     \item{ks_test}{List with \code{statistic} and \code{p_value}
#'       from the Kolmogorov-Smirnov test for uniformity.}
#'     \item{n}{Integer; number of forecast-observation pairs used.}
#'   }
#'
#' @details
#' The PIT value for a single forecast-observation pair is defined as
#' the quantile of the observation within the forecast distribution.
#' Under the assumption of Gaussian prediction intervals (which the
#' parametric simulation in \code{\link{forecast}} approximates),
#' the PIT is computed as
#' \eqn{\Phi((y - \hat{y}) / \hat{\sigma})}
#' where \eqn{\hat{y}} is the predicted median, \eqn{\hat{\sigma}}
#' is derived from the prediction interval width, and \eqn{\Phi} is
#' the standard normal CDF.
#'
#' The reliability diagram plots observed coverage against nominal
#' coverage at levels 10 through 90 percent. Perfect calibration
#' lies on the diagonal.
#'
#' @references
#' Gneiting T, Balabdaoui F, Raftery AE (2007). Probabilistic
#' forecasts, calibration and sharpness. \emph{Journal of the Royal
#' Statistical Society: Series B}, 69(2), 243--268.
#' \doi{10.1111/j.1467-9868.2007.00587.x}
#'
#' @seealso \code{\link{recalibrate}} for post-hoc recalibration,
#'   \code{\link{score_forecasts}} for proper scoring rules.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 20, seed = 1)
#' bt <- backtest(sim, engines = "mlr",
#'   horizons = c(7, 14), min_train = 42)
#' cal <- calibrate(bt)
#' cal
#' }
#'
#' @export
calibrate <- function(x, observed = NULL, n_bins = 10L) {
  UseMethod("calibrate")
}


#' @rdname calibrate
#' @export
calibrate.lfq_backtest <- function(x, observed = NULL, n_bins = 10L) {

  df <- x[!is.na(x$observed) & !is.na(x$predicted), ]

  if (nrow(df) == 0L) {
    cli::cli_abort("No valid forecast-observation pairs for calibration.")
  }

  pit <- .compute_pit(df$observed, df$predicted, df$lower, df$upper)
  .build_calibration_report(pit, n_bins)
}


#' @rdname calibrate
#' @export
calibrate.lfq_forecast <- function(x, observed = NULL, n_bins = 10L) {

  fc_rows <- x[x$.type == "forecast", ]

  if (is.null(observed)) {
    cli::cli_abort(
      "Supply {.arg observed} frequencies for {.cls lfq_forecast} objects."
    )
  }

  if (length(observed) != nrow(fc_rows)) {
    cli::cli_abort(
      "{.arg observed} length ({length(observed)}) must match forecast rows ({nrow(fc_rows)})."
    )
  }

  pit <- .compute_pit(observed, fc_rows$.median, fc_rows$.lower,
                      fc_rows$.upper)
  .build_calibration_report(pit, n_bins)
}


#' Compute PIT values from prediction intervals
#'
#' Assumes the forecast distribution is approximately Gaussian,
#' centred on the predicted median with scale derived from the
#' prediction interval width.
#' @noRd
.compute_pit <- function(observed, predicted, lower, upper) {
  # Derive implied standard deviation from the 95% interval
  # width = 2 * z_{0.975} * sigma => sigma = width / (2 * 1.96)
  width <- upper - lower
  sigma <- pmax(width / (2 * stats::qnorm(0.975)), 1e-10)
  z <- (observed - predicted) / sigma
  stats::pnorm(z)
}


#' Build calibration report from PIT values
#' @noRd
.build_calibration_report <- function(pit, n_bins) {

  n <- length(pit)

  # PIT histogram
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  counts <- as.integer(table(cut(pit, breaks = breaks,
                                  include.lowest = TRUE)))
  pit_hist <- tibble::tibble(
    bin      = seq_len(n_bins),
    count    = counts,
    density  = counts / n,
    expected = 1 / n_bins
  )

  # Reliability diagram: observed coverage at nominal levels
  nominal_levels <- seq(0.1, 0.9, by = 0.1)
  observed_cov <- vapply(nominal_levels, function(level) {
    # What fraction of PIT values fall in the central 'level' interval?
    alpha <- (1 - level) / 2
    mean(pit >= alpha & pit <= 1 - alpha)
  }, numeric(1L))

  reliability <- tibble::tibble(
    nominal  = nominal_levels,
    observed = observed_cov
  )

  # KS test for uniformity
  ks <- stats::ks.test(pit, "punif")

  structure(
    list(
      pit_values    = pit,
      pit_histogram = pit_hist,
      reliability   = reliability,
      ks_test       = list(statistic = unname(ks$statistic),
                           p_value   = ks$p.value),
      n             = n
    ),
    class = "calibration_report"
  )
}


#' @export
print.calibration_report <- function(x, ...) {
  cli::cli_h3("Calibration report")
  cli::cli_text("{x$n} forecast-observation pairs")
  cli::cli_text(
    "KS test for PIT uniformity: D = {round(x$ks_test$statistic, 4)}, p = {format.pval(x$ks_test$p_value, digits = 3)}"
  )

  # Summary of coverage calibration
  rel <- x$reliability
  mean_err <- mean(abs(rel$observed - rel$nominal))
  cli::cli_text(
    "Mean absolute calibration error: {round(mean_err, 4)}"
  )

  invisible(x)
}


#' Plot calibration diagnostics
#'
#' Produces either a reliability diagram (default) or a PIT histogram.
#'
#' @param x A \code{calibration_report} object.
#' @param type Which panel to display: \code{"reliability"} (default)
#'   or \code{"pit"} for the PIT histogram.
#' @param ... Unused.
#'
#' @return A ggplot object.
#'
#' @export


#' @rdname plot.calibration_report
#' @param type Which panel to display: \code{"reliability"} (default)
#'   or \code{"pit"} for the PIT histogram.
#' @export
plot.calibration_report <- function(x, type = c("reliability", "pit"),
                                    ...) {

  type <- match.arg(type)

  if (type == "pit") {
    ggplot2::ggplot(x$pit_histogram,
      ggplot2::aes(x = .data$bin, y = .data$density)) +
      ggplot2::geom_col(fill = "#4292c6", alpha = 0.8) +
      ggplot2::geom_hline(yintercept = x$pit_histogram$expected[1],
                          linetype = "dashed", colour = "grey40") +
      ggplot2::labs(x = "PIT bin", y = "Density",
                   title = "PIT histogram") +
      ggplot2::theme_minimal(base_size = 11)
  } else {
    ggplot2::ggplot(x$reliability,
      ggplot2::aes(x = .data$nominal, y = .data$observed)) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", colour = "grey40") +
      ggplot2::geom_point(size = 2.5, colour = "#d6604d") +
      ggplot2::geom_line(colour = "#d6604d") +
      ggplot2::scale_x_continuous(limits = c(0, 1)) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = "Nominal coverage", y = "Observed coverage",
                   title = "Reliability diagram") +
      ggplot2::theme_minimal(base_size = 11)
  }
}
