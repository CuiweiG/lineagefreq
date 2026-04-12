#' Score backtest forecast accuracy
#'
#' Computes standardized accuracy metrics from backtesting results.
#'
#' @param bt An `lfq_backtest` object from [backtest()].
#' @param metrics Character vector of metrics to compute:
#'   * `"mae"`: Mean absolute error of frequency.
#'   * `"rmse"`: Root mean squared error.
#'   * `"coverage"`: Proportion within prediction intervals.
#'   * `"wis"`: Simplified weighted interval score for the single
#'     prediction interval stored in the backtest (typically 95%).
#'   * `"crps"`: Continuous Ranked Probability Score, assuming
#'     Gaussian forecast distribution (Gneiting and Raftery, 2007).
#'   * `"log_score"`: Logarithmic scoring rule evaluated at the
#'     observed value under the Gaussian forecast density.
#'   * `"dss"`: Dawid-Sebastiani Score, a proper scoring rule based
#'     on the predictive mean and variance.
#'   * `"calibration"`: Mean squared calibration error across
#'     nominal coverage levels 10\%--90\%.
#'
#' @return A tibble with columns: `engine`, `horizon`, `metric`,
#'   `value`.
#'
#' @references
#' Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating
#' epidemic forecasts in an interval format. \emph{PLoS
#' Computational Biology}, 17(2):e1008618.
#' \doi{10.1371/journal.pcbi.1008618}
#'
#' Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
#' prediction, and estimation. \emph{Journal of the American
#' Statistical Association}, 102(477), 359--378.
#' \doi{10.1198/016214506000001437}
#'
#' @seealso [compare_models()] to rank engines based on scores.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 20, seed = 1)
#' bt <- backtest(sim, engines = "mlr",
#'   horizons = c(7, 14), min_train = 42)
#' score_forecasts(bt)
#' }
#'
#' @export
score_forecasts <- function(bt,
                            metrics = c("mae", "rmse", "coverage",
                                        "wis", "crps", "log_score",
                                        "dss", "calibration")) {

  if (!inherits(bt, "lfq_backtest")) {
    cli::cli_abort("{.arg bt} must be an {.cls lfq_backtest} object.")
  }

  metrics <- match.arg(metrics, several.ok = TRUE)

  if (nrow(bt) == 0L) {
    cli::cli_warn("No backtest results to score.")
    return(tibble::tibble(engine = character(), horizon = integer(),
                          metric = character(), value = numeric()))
  }

  bt <- bt[!is.na(bt$observed) & !is.na(bt$predicted), ]

  score_rows <- list()

  for (eng in unique(bt$engine)) {
    for (h in unique(bt$horizon)) {
      sub <- bt[bt$engine == eng & bt$horizon == h, ]
      if (nrow(sub) == 0L) next

      for (m in metrics) {
        val <- switch(m,
          mae         = mean(abs(sub$predicted - sub$observed)),
          rmse        = sqrt(mean((sub$predicted - sub$observed)^2)),
          coverage    = mean(sub$observed >= sub$lower &
                               sub$observed <= sub$upper, na.rm = TRUE),
          wis         = .compute_wis(sub),
          crps        = .compute_crps(sub),
          log_score   = .compute_log_score(sub),
          dss         = .compute_dss(sub),
          calibration = .compute_calibration_score(sub)
        )
        score_rows <- c(score_rows, list(tibble::tibble(
          engine  = eng,
          horizon = h,
          metric  = m,
          value   = val
        )))
      }
    }
  }

  dplyr::bind_rows(score_rows)
}


#' Simplified Weighted Interval Score (internal)
#'
#' Computes WIS for a single interval level (95% by default).
#' The full WIS (Bracher et al. 2021) averages over multiple
#' quantile levels; this simplified version uses only the
#' prediction interval stored in the backtest output.
#' @noRd
.compute_wis <- function(df) {
  alpha      <- 0.05
  overshoot  <- pmax(df$lower - df$observed, 0)
  undershoot <- pmax(df$observed - df$upper, 0)
  width      <- df$upper - df$lower
  wis        <- (alpha / 2) * width + overshoot + undershoot
  mean(wis, na.rm = TRUE)
}


#' Gaussian CRPS (internal)
#'
#' Continuous Ranked Probability Score under a Gaussian forecast
#' distribution with mean = predicted and sigma derived from the
#' prediction interval width. Closed-form: Gneiting & Raftery 2007.
#' @noRd
.compute_crps <- function(df) {
  sigma <- .implied_sigma(df)
  z <- (df$observed - df$predicted) / sigma
  # CRPS = sigma * (z * (2*Phi(z) - 1) + 2*phi(z) - 1/sqrt(pi))
  crps <- sigma * (z * (2 * stats::pnorm(z) - 1) +
                     2 * stats::dnorm(z) - 1 / sqrt(pi))
  mean(crps, na.rm = TRUE)
}


#' Gaussian log score (internal)
#'
#' Negative log-likelihood of the observation under a Gaussian
#' forecast distribution. Lower is better.
#' @noRd
.compute_log_score <- function(df) {
  sigma <- .implied_sigma(df)
  # -log(dnorm(observed, predicted, sigma))
  nll <- 0.5 * log(2 * pi) + log(sigma) +
    0.5 * ((df$observed - df$predicted) / sigma)^2
  mean(nll, na.rm = TRUE)
}


#' Dawid-Sebastiani Score (internal)
#'
#' DSS = log(sigma^2) + ((y - mu) / sigma)^2. A proper scoring
#' rule that depends only on the first two moments of the forecast
#' distribution.
#' @noRd
.compute_dss <- function(df) {
  sigma <- .implied_sigma(df)
  dss <- log(sigma^2) + ((df$observed - df$predicted) / sigma)^2
  mean(dss, na.rm = TRUE)
}


#' Mean squared calibration error (internal)
#'
#' Computes coverage at nominal levels 10%-90% and returns the
#' mean squared difference from nominal.
#' @noRd
.compute_calibration_score <- function(df) {
  sigma <- .implied_sigma(df)
  z <- (df$observed - df$predicted) / sigma
  nominal <- seq(0.1, 0.9, by = 0.1)
  msce <- mean(vapply(nominal, function(level) {
    z_crit <- stats::qnorm(1 - (1 - level) / 2)
    obs_cov <- mean(abs(z) <= z_crit, na.rm = TRUE)
    (obs_cov - level)^2
  }, numeric(1L)))
  msce
}


#' Implied sigma from prediction interval (internal)
#'
#' Floors sigma at 0.001 to prevent numerical overflow in log score
#' and DSS when prediction intervals have near-zero width (e.g.,
#' fitted values with no uncertainty in the backtest).
#' @noRd
.implied_sigma <- function(df) {
  width <- df$upper - df$lower
  pmax(width / (2 * stats::qnorm(0.975)), 0.001)
}


#' Compare model engines from backtest scores
#'
#' Summarises and ranks engines across horizons based on forecast
#' accuracy scores.
#'
#' @param scores Output of [score_forecasts()].
#' @param by Grouping variable(s). Default `"engine"`.
#'
#' @return A tibble with average scores per group, sorted by MAE.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 20, seed = 1)
#' bt <- backtest(sim, engines = "mlr",
#'   horizons = c(7, 14), min_train = 42)
#' sc <- score_forecasts(bt)
#' compare_models(sc)
#' }
#'
#' @export
compare_models <- function(scores, by = "engine") {

  if (!is.data.frame(scores) || nrow(scores) == 0L) {
    cli::cli_warn("No scores to compare.")
    return(tibble::tibble())
  }

  wide <- scores |>
    tidyr::pivot_wider(
      id_cols    = dplyr::all_of(by),
      names_from = "metric",
      values_from = "value",
      values_fn  = mean
    )

  sort_col <- if ("mae" %in% names(wide)) "mae" else {
    num_cols <- names(wide)[vapply(wide, is.numeric, logical(1L))]
    if (length(num_cols) > 0L) num_cols[1L] else NULL
  }

  if (!is.null(sort_col)) {
    wide <- dplyr::arrange(wide, .data[[sort_col]])
  }

  wide
}


#' Plot backtest scores
#'
#' Creates a panel plot of forecast accuracy by engine and horizon.
#'
#' @param scores Output of [score_forecasts()].
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 20, seed = 1)
#' bt <- backtest(sim, engines = "mlr",
#'   horizons = c(7, 14), min_train = 42)
#' sc <- score_forecasts(bt)
#' plot_backtest(sc)
#' }
#'
#' @export
plot_backtest <- function(scores) {

  if (!is.data.frame(scores) || nrow(scores) == 0L) {
    cli::cli_abort("No scores to plot.")
  }

  ggplot2::ggplot(scores, ggplot2::aes(
    x    = factor(.data$horizon),
    y    = .data$value,
    fill = .data$engine
  )) +
    ggplot2::geom_col(position = "dodge", alpha = 0.8) +
    ggplot2::facet_wrap(~ .data$metric, scales = "free_y") +
    ggplot2::labs(
      x    = "Forecast horizon (days)",
      y    = "Score",
      fill = "Engine"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
