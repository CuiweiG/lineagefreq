#' Score backtest forecast accuracy
#'
#' Computes standardized accuracy metrics from backtesting results.
#'
#' @param bt An [lfq_backtest] object from [backtest()].
#' @param metrics Character vector of metrics to compute:
#'   * `"mae"`: Mean absolute error of frequency.
#'   * `"rmse"`: Root mean squared error.
#'   * `"coverage"`: Proportion within prediction intervals.
#'   * `"wis"`: Weighted interval score.
#'
#' @return A tibble with columns: `engine`, `horizon`, `metric`,
#'   `value`.
#'
#' @seealso [compare_models()] to rank engines based on scores.
#'
#' @examples
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 20, seed = 1)
#' bt <- backtest(sim, engines = "mlr",
#'   horizons = c(7, 14), min_train = 42)
#' score_forecasts(bt)
#'
#' @export
score_forecasts <- function(bt,
                            metrics = c("mae", "rmse", "coverage",
                                        "wis")) {

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
          mae      = mean(abs(sub$predicted - sub$observed)),
          rmse     = sqrt(mean((sub$predicted - sub$observed)^2)),
          coverage = mean(sub$observed >= sub$lower &
                            sub$observed <= sub$upper, na.rm = TRUE),
          wis      = .compute_wis(sub)
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


#' Weighted Interval Score (internal)
#' @noRd
.compute_wis <- function(df) {
  alpha      <- 0.05
  overshoot  <- pmax(df$lower - df$observed, 0)
  undershoot <- pmax(df$observed - df$upper, 0)
  width      <- df$upper - df$lower
  wis        <- (alpha / 2) * width + overshoot + undershoot
  mean(wis, na.rm = TRUE)
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
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 20, seed = 1)
#' bt <- backtest(sim, engines = "mlr",
#'   horizons = c(7, 14), min_train = 42)
#' sc <- score_forecasts(bt)
#' compare_models(sc)
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
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 20, seed = 1)
#' bt <- backtest(sim, engines = "mlr",
#'   horizons = c(7, 14), min_train = 42)
#' sc <- score_forecasts(bt)
#' plot_backtest(sc)
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
