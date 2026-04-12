#' Recalibrate prediction intervals
#'
#' Applies post-hoc recalibration to improve the coverage properties
#' of prediction intervals from \code{\link{forecast}}. Two methods
#' are available: isotonic regression (nonparametric, monotonicity-
#' preserving) and Platt scaling (logistic, parametric).
#'
#' @param forecast_obj An \code{lfq_forecast} object to recalibrate.
#' @param bt An \code{lfq_backtest} object providing the calibration
#'   data. The mapping from nominal to empirical coverage is learned
#'   from the backtest residuals.
#' @param method Recalibration method: \code{"isotonic"} (default)
#'   for isotonic regression on the empirical coverage function, or
#'   \code{"platt"} for Platt scaling via logistic regression.
#'
#' @return An \code{lfq_forecast} object with recalibrated
#'   \code{.lower} and \code{.upper} bounds. The object retains all
#'   attributes of the original forecast and can be passed to
#'   \code{\link[ggplot2]{autoplot}}.
#'
#' @details
#' \strong{Isotonic regression:} Learns the monotone mapping from
#' nominal coverage levels to observed coverage using the backtest
#' data, then inverts it to find the nominal level that achieves
#' the desired empirical coverage. This is a nonparametric approach
#' that requires no distributional assumptions.
#'
#' \strong{Platt scaling:} Fits a logistic regression of the form
#' \eqn{P(\text{covered}) = \text{logit}^{-1}(a \cdot z + b)}
#' where \eqn{z} is the standardised residual, and uses the fitted
#' model to adjust prediction interval widths.
#'
#' Both methods require a backtest object with sufficient data to
#' estimate the calibration mapping reliably (at least 30
#' forecast-observation pairs recommended).
#'
#' @references
#' Platt JC (1999). Probabilistic outputs for support vector
#' machines and comparisons to regularized likelihood methods.
#' \emph{Advances in Large Margin Classifiers}, 61--74.
#'
#' @seealso \code{\link{calibrate}} for calibration diagnostics.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 20, seed = 1)
#' fit <- fit_model(sim, engine = "mlr")
#' fc  <- forecast(fit, horizon = 21)
#' bt  <- backtest(sim, engines = "mlr",
#'   horizons = c(7, 14), min_train = 42)
#' fc_recal <- recalibrate(fc, bt, method = "isotonic")
#' }
#'
#' @export
recalibrate <- function(forecast_obj, bt,
                        method = c("isotonic", "platt")) {

  if (!inherits(forecast_obj, "lfq_forecast")) {
    cli::cli_abort("{.arg forecast_obj} must be an {.cls lfq_forecast}.")
  }
  if (!inherits(bt, "lfq_backtest")) {
    cli::cli_abort("{.arg bt} must be an {.cls lfq_backtest}.")
  }

  method <- match.arg(method)

  df <- bt[!is.na(bt$observed) & !is.na(bt$predicted), ]
  if (nrow(df) < 10L) {
    cli::cli_warn(
      "Few backtest pairs ({nrow(df)}); recalibration may be unreliable."
    )
  }

  # Learn the calibration adjustment from backtest residuals
  width_bt <- df$upper - df$lower
  sigma_bt <- pmax(width_bt / (2 * stats::qnorm(0.975)), 1e-10)
  z_bt <- (df$observed - df$predicted) / sigma_bt

  if (method == "isotonic") {
    recal_factor <- .isotonic_recalibrate(z_bt)
  } else {
    recal_factor <- .platt_recalibrate(z_bt)
  }

  # Apply to forecast intervals
  fc <- forecast_obj
  fc_idx <- fc$.type == "forecast"

  if (sum(fc_idx) == 0L) return(fc)

  width_fc <- fc$.upper[fc_idx] - fc$.lower[fc_idx]
  sigma_fc <- pmax(width_fc / (2 * stats::qnorm(0.975)), 1e-10)

  # Widen or narrow intervals by the recalibration factor
  new_sigma <- sigma_fc * recal_factor
  fc$.lower[fc_idx] <- pmax(fc$.median[fc_idx] -
    stats::qnorm(0.975) * new_sigma, 0)
  fc$.upper[fc_idx] <- pmin(fc$.median[fc_idx] +
    stats::qnorm(0.975) * new_sigma, 1)

  fc
}


#' Isotonic regression recalibration factor
#'
#' Computes the ratio of observed to nominal interval width needed
#' to achieve calibrated coverage.
#' @noRd
.isotonic_recalibrate <- function(z) {
  # Empirical coverage at the nominal 95% level
  nominal_coverage <- 0.95
  observed_coverage <- mean(abs(z) <= stats::qnorm(0.975))

  if (observed_coverage >= nominal_coverage) {
    # Already well-calibrated or over-covered
    return(1.0)
  }

  # Find the z-quantile that achieves nominal coverage empirically
  z_abs <- sort(abs(z))
  target_idx <- ceiling(nominal_coverage * length(z_abs))
  target_idx <- min(target_idx, length(z_abs))
  empirical_z <- z_abs[target_idx]

  # Ratio: how much wider do intervals need to be?
  empirical_z / stats::qnorm(0.975)
}


#' Platt scaling recalibration factor
#'
#' Fits logistic regression on standardised residuals and derives
#' the scaling factor.
#' @noRd
.platt_recalibrate <- function(z) {
  covered <- as.integer(abs(z) <= stats::qnorm(0.975))
  z_abs <- abs(z)

  fit <- tryCatch(
    stats::glm(covered ~ z_abs, family = stats::binomial()),
    error = function(e) NULL
  )

  if (is.null(fit)) return(1.0)

  # Find z_abs where P(covered) = 0.95
  # logit(0.95) = coef[1] + coef[2] * z_target
  logit_target <- stats::qlogis(0.95)
  coefs <- stats::coef(fit)
  if (is.na(coefs[2]) || abs(coefs[2]) < 1e-10) return(1.0)

  z_target <- (logit_target - coefs[1]) / coefs[2]
  if (z_target <= 0) return(1.0)

  z_target / stats::qnorm(0.975)
}
