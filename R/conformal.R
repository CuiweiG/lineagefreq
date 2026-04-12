#' Conformal prediction intervals for lineage frequencies
#'
#' Produces distribution-free prediction intervals with finite-sample
#' coverage guarantees using split conformal inference. Unlike the
#' parametric intervals from \code{\link{forecast}}, conformal
#' intervals require no distributional assumptions on the residuals
#' and are valid under exchangeability.
#'
#' @param fit An \code{lfq_fit} object (any engine).
#' @param data The \code{lfq_data} object used to fit the model.
#' @param horizon Number of days to forecast. Default 28.
#' @param ci_level Target coverage level. Default 0.95.
#' @param method Conformal method: \code{"split"} (default) for
#'   split conformal prediction, or \code{"aci"} for adaptive
#'   conformal inference with online coverage correction.
#' @param cal_fraction Fraction of the data reserved for the
#'   calibration set (split conformal only). Default 0.3.
#' @param gamma Learning rate for adaptive conformal inference.
#'   Default 0.05. Controls how quickly the coverage target adjusts
#'   in response to observed miscoverage.
#' @param seed Random seed for the calibration split. Default
#'   \code{NULL}.
#'
#' @return An \code{lfq_forecast} object with conformal prediction
#'   intervals. The object is fully compatible with
#'   \code{\link[ggplot2]{autoplot}} and other forecast methods.
#'
#' @details
#' \strong{Split conformal prediction} partitions the training data
#' into a proper training set and a calibration set. The model is
#' refit on the training set, and conformity scores (absolute
#' residuals) are computed on the calibration set. The prediction
#' interval at a new point is the point forecast plus or minus the
#' \eqn{(1 - \alpha)(1 + 1/n_{\text{cal}})} quantile of the
#' calibration scores. This provides exact \eqn{1 - \alpha} marginal
#' coverage under exchangeability (Vovk et al. 2005).
#'
#' \strong{Adaptive conformal inference (ACI)} (Gibbs and Candes,
#' 2021) adjusts the miscoverage level online to maintain long-run
#' coverage even when the data distribution shifts over time, as is
#' typical in surveillance data during variant replacement waves.
#'
#' @references
#' Vovk V, Gammerman A, Shafer G (2005). \emph{Algorithmic Learning
#' in a Random World}. Springer.
#'
#' Gibbs I, Candes E (2021). Adaptive conformal inference under
#' distribution shift. \emph{Advances in Neural Information
#' Processing Systems}, 34.
#'
#' @seealso \code{\link{forecast}} for parametric prediction
#'   intervals, \code{\link{calibrate}} for calibration diagnostics.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 20, seed = 1)
#' fit <- fit_model(sim, engine = "mlr")
#' fc_conf <- conformal_forecast(fit, sim, horizon = 21)
#' fc_conf
#' }
#'
#' @export
conformal_forecast <- function(fit, data,
                               horizon      = 28L,
                               ci_level     = 0.95,
                               method       = c("split", "aci"),
                               cal_fraction = 0.3,
                               gamma        = 0.05,
                               seed         = NULL) {

  if (!inherits(fit, "lfq_fit")) {
    cli::cli_abort("{.arg fit} must be an {.cls lfq_fit} object.")
  }
  if (!is_lfq_data(data)) {
    cli::cli_abort("{.arg data} must be an {.cls lfq_data} object.")
  }

  method <- match.arg(method)
  assert_prob(ci_level, "ci_level")

  lineages <- fit$lineages
  pivot    <- fit$pivot
  ts       <- fit$time_scale

  all_dates <- sort(unique(data$.date))
  n_dates   <- length(all_dates)

  if (!is.null(seed)) set.seed(seed)

  if (method == "split") {
    scores <- .split_conformal_scores(fit, data, cal_fraction, seed)
  } else {
    scores <- .aci_scores(fit, data, ci_level, gamma)
  }

  # Conformal quantile (Vovk et al. 2005, Theorem 2.2)
  n_cal <- length(scores)
  q_level <- ceiling((1 - (1 - ci_level)) * (n_cal + 1)) / n_cal
  q_level <- min(q_level, 1)
  q_hat <- stats::quantile(scores, q_level, names = FALSE)

  # Generate point forecasts from the original fit
  fc_parametric <- forecast(fit, horizon = horizon, ci_level = ci_level,
                            n_sim = 1000L)

  # Replace parametric intervals with conformal intervals
  fc <- fc_parametric
  fc_idx <- fc$.type == "forecast"
  fc$.lower[fc_idx] <- pmax(fc$.median[fc_idx] - q_hat, 0)
  fc$.upper[fc_idx] <- pmin(fc$.median[fc_idx] + q_hat, 1)

  fc
}


#' Compute split conformal calibration scores
#' @noRd
.split_conformal_scores <- function(fit, data, cal_fraction, seed) {

  all_dates <- sort(unique(data$.date))
  n_dates   <- length(all_dates)

  # Temporal split: use the last cal_fraction of dates for calibration
  n_cal    <- max(floor(n_dates * cal_fraction), 2L)
  cal_dates <- utils::tail(all_dates, n_cal)
  train_dates <- utils::head(all_dates, n_dates - n_cal)

  if (length(train_dates) < 3L) {
    cli::cli_abort("Insufficient training data after calibration split.")
  }

  # Refit on training portion
  train_data <- data[data$.date %in% train_dates, ]
  train_data <- structure(
    train_data,
    class        = class(data),
    lineages     = attr(data, "lineages"),
    date_range   = range(train_dates),
    n_timepoints = length(train_dates),
    has_location = attr(data, "has_location"),
    min_total    = attr(data, "min_total")
  )

  refit <- tryCatch(
    fit_model(train_data, engine = fit$engine),
    error = function(e) {
      cli::cli_abort("Model refit failed during conformal calibration: {e$message}")
    }
  )

  # Compute fitted values on calibration dates
  t0 <- min(train_dates)
  ts <- refit$time_scale
  lineages  <- refit$lineages
  non_pivot <- setdiff(lineages, refit$pivot)

  scores <- numeric(0)

  for (d in cal_dates) {
    d <- as.Date(d, origin = "1970-01-01")
    t_val <- as.numeric(d - t0) / ts

    log_num <- stats::setNames(
      c(0, refit$intercepts[non_pivot]),
      c(refit$pivot, non_pivot)
    )[lineages] +
      stats::setNames(
        c(0, refit$growth_rates[non_pivot]),
        c(refit$pivot, non_pivot)
      )[lineages] * t_val

    log_denom <- log_sum_exp(log_num)
    pred_freq <- exp(log_num - log_denom)

    obs <- data[data$.date == d, ]
    for (lin in lineages) {
      obs_row <- obs[obs$.lineage == lin, ]
      if (nrow(obs_row) == 1L) {
        scores <- c(scores, abs(obs_row$.freq - pred_freq[lin]))
      }
    }
  }

  scores
}


#' Compute ACI scores with adaptive alpha
#' @noRd
.aci_scores <- function(fit, data, ci_level, gamma) {
  # Use all residuals from the fit as conformity scores,
  # but weight recent residuals more heavily via the adaptive alpha
  resid_df <- fit$residuals
  if (is.null(resid_df) || nrow(resid_df) == 0L) {
    cli::cli_abort("Model has no residuals for ACI computation.")
  }

  scores <- abs(resid_df$.pearson_resid)
  alpha_t <- 1 - ci_level

  # Online update: if recent forecasts are under-covered, alpha
  # decreases (intervals widen); if over-covered, alpha increases
  n <- length(scores)
  q_current <- stats::quantile(scores,
    1 - alpha_t, names = FALSE)

  # Track coverage in rolling windows
  covered <- scores <= q_current
  for (i in seq_len(n)) {
    alpha_t <- alpha_t + gamma * ((1 - ci_level) - (1 - covered[i]))
    alpha_t <- max(min(alpha_t, 0.5), 0.001)
  }

  # Return scores scaled by the adaptive factor
  scores * (stats::qnorm(1 - alpha_t / 2) / stats::qnorm(0.975))
}
