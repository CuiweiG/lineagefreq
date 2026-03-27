#' Forecast lineage frequencies (generic)
#'
#' @param object A model object.
#' @param ... Additional arguments passed to methods.
#' @return A forecast object.
#' @export
forecast <- function(object, ...) UseMethod("forecast")


#' Forecast lineage frequencies
#'
#' Projects lineage frequencies forward in time using the fitted model.
#' Prediction uncertainty is quantified by parametric simulation from
#' the estimated parameter distribution.
#'
#' @param object An `lfq_fit` object.
#' @param horizon Number of days to forecast. Default 28 (4 weeks).
#' @param ci_level Confidence level for prediction intervals.
#'   Default 0.95.
#' @param n_sim Number of parameter draws for prediction intervals.
#'   Default 1000.
#' @param ... Unused.
#'
#' @return An `lfq_forecast` object (tibble subclass) with columns:
#' \describe{
#'   \item{.date}{Date.}
#'   \item{.lineage}{Lineage name.}
#'   \item{.median}{Median predicted frequency.}
#'   \item{.lower}{Lower prediction bound.}
#'   \item{.upper}{Upper prediction bound.}
#'   \item{.type}{"fitted" or "forecast".}
#' }
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
#' fit <- fit_model(sim, engine = "mlr")
#' fc <- forecast(fit, horizon = 21)
#' fc
#' }
#'
#' @export
forecast.lfq_fit <- function(object, horizon  = 28L,
                             ci_level = 0.95,
                             n_sim    = 1000L, ...) {

  assert_pos_int(horizon, "horizon")
  assert_prob(ci_level, "ci_level")
  assert_pos_int(n_sim, "n_sim")

  lineages  <- object$lineages
  pivot     <- object$pivot
  non_pivot <- setdiff(lineages, pivot)
  n_lin     <- length(lineages)
  ts        <- object$time_scale

  # --- Fitted part: convert fitted_values to standard format ---
  fitted_df <- object$fitted_values
  fitted_df$.lower  <- fitted_df$.fitted_freq
  fitted_df$.upper  <- fitted_df$.fitted_freq
  fitted_df         <- dplyr::rename(fitted_df, .median = ".fitted_freq")
  fitted_df$.type   <- "fitted"

  # --- Future dates at same interval as original data ---
  max_date <- max(object$date_range)
  t0       <- min(object$date_range)

  dates_sorted      <- sort(unique(object$fitted_values$.date))
  original_interval <- if (length(dates_sorted) > 1L) {
    as.integer(stats::median(diff(as.numeric(dates_sorted))))
  } else {
    7L
  }

  future_dates <- seq(max_date + original_interval,
                      max_date + horizon,
                      by = original_interval)
  if (length(future_dates) == 0L) {
    future_dates <- max_date + as.integer(horizon)
  }

  # --- Parameter simulation via MVN ---
  par_est <- c(
    unname(object$intercepts[non_pivot]),
    unname(object$growth_rates[non_pivot])
  )
  vcov_mat <- object$vcov_matrix

  # Ensure vcov is positive-definite
  vcov_safe <- tryCatch({
    chol(vcov_mat)   # test; throws if not PD
    vcov_mat
  }, error = function(e) {
    ev <- eigen(vcov_mat, symmetric = TRUE)
    ev$values <- pmax(ev$values, 1e-8)
    ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
  })

  par_draws <- MASS::mvrnorm(n_sim, mu = par_est, Sigma = vcov_safe)

  alpha_q <- (1 - ci_level) / 2

  # --- Forecast rows ---
  forecast_rows <- vector("list", length(future_dates) * n_lin)
  idx <- 0L

  for (j in seq_along(future_dates)) {
    t_val <- as.numeric(future_dates[j] - t0) / ts

    # n_sim × n_lin frequency matrix
    freq_matrix <- matrix(NA_real_, nrow = n_sim, ncol = n_lin)

    for (s in seq_len(n_sim)) {
      alpha_draw <- par_draws[s, seq_len(n_lin - 1L)]
      delta_draw <- par_draws[s, (n_lin - 1L) + seq_len(n_lin - 1L)]

      alpha_full <- stats::setNames(
        c(0, alpha_draw), c(pivot, non_pivot)
      )[lineages]
      delta_full <- stats::setNames(
        c(0, delta_draw), c(pivot, non_pivot)
      )[lineages]

      log_num   <- alpha_full + delta_full * t_val
      log_denom <- log_sum_exp(log_num)
      freq_matrix[s, ] <- exp(log_num - log_denom)
    }

    for (v_idx in seq_along(lineages)) {
      vals <- freq_matrix[, v_idx]
      idx  <- idx + 1L
      forecast_rows[[idx]] <- tibble::tibble(
        .date    = future_dates[j],
        .lineage = lineages[v_idx],
        .median  = stats::median(vals),
        .lower   = stats::quantile(vals, alpha_q,     names = FALSE),
        .upper   = stats::quantile(vals, 1 - alpha_q, names = FALSE),
        .type    = "forecast"
      )
    }
  }

  forecast_df <- dplyr::bind_rows(forecast_rows)

  # --- Combine fitted + forecast ---
  combined <- dplyr::bind_rows(fitted_df, forecast_df)
  combined <- dplyr::arrange(combined, .data$.date, .data$.lineage)

  structure(
    combined,
    class          = c("lfq_forecast", class(tibble::tibble())),
    horizon        = as.integer(horizon),
    ci_level       = ci_level,
    engine         = object$engine,
    pivot          = pivot,
    lineages       = lineages,
    forecast_start = max_date + 1L
  )
}


#' @return The input object, invisibly.
#' @export
print.lfq_forecast <- function(x, ...) {
  fc_start <- attr(x, "forecast_start")
  horizon  <- attr(x, "horizon")
  ci       <- attr(x, "ci_level")
  n_fc     <- sum(x$.type == "forecast")
  n_fit    <- sum(x$.type == "fitted")

  cli::cli_h3("Lineage frequency forecast")
  cli::cli_text("Engine: {attr(x, 'engine')}")
  cli::cli_text(
    "Forecast start: {fc_start}  |  Horizon: {horizon} days"
  )
  cli::cli_text("CI level: {ci * 100}%")
  cli::cli_text("{n_fit} fitted + {n_fc} forecast rows")
  cat("\n")
  NextMethod()
}
