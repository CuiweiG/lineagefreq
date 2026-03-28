#' Pipe-friendly model fitting
#'
#' Enables tidyverse-style chaining:
#' `data |> lfq_fit("mlr") |> lfq_forecast(28) |> lfq_score()`
#'
#' @param data An [lfq_data] object.
#' @param engine Engine name. Default `"mlr"`.
#' @param ... Passed to [fit_model()].
#'
#' @return An `lfq_fit` object.
#'
#' @examples
#' data(sarscov2_us_2022)
#' sarscov2_us_2022 |>
#'   lfq_data(lineage = variant, date = date, count = count, total = total) |>
#'   lfq_fit("mlr") |>
#'   lfq_advantage(generation_time = 5)
#'
#' @export
lfq_fit <- function(data, engine = "mlr", ...) {
  fit_model(data, engine = engine, ...)
}

#' Pipe-friendly growth advantage extraction
#'
#' @param fit An `lfq_fit` object.
#' @param type Output type. Default `"relative_Rt"`.
#' @param generation_time Mean generation time in days.
#' @param ... Passed to [growth_advantage()].
#'
#' @return A tibble of growth advantages.
#'
#' @export
lfq_advantage <- function(fit, type = "relative_Rt",
                          generation_time = NULL, ...) {
  growth_advantage(fit, type = type,
                   generation_time = generation_time, ...)
}

#' Pipe-friendly forecasting
#'
#' @param fit An `lfq_fit` object.
#' @param horizon Forecast horizon in days. Default 28.
#' @param ... Passed to [forecast()].
#'
#' @return An `lfq_forecast` object.
#'
#' @export
lfq_forecast <- function(fit, horizon = 28L, ...) {
  forecast(fit, horizon = horizon, ...)
}

#' Pipe-friendly backtesting + scoring
#'
#' Runs backtest and returns scores in one step.
#'
#' @param data An [lfq_data] object.
#' @param engines Character vector of engine names.
#' @param horizons Forecast horizons in days.
#' @param metrics Score metrics.
#' @param ... Passed to [backtest()].
#'
#' @return A tibble of scores.
#'
#' @export
lfq_score <- function(data, engines = "mlr",
                      horizons = c(14, 28),
                      metrics = c("mae", "coverage"),
                      ...) {
  bt <- backtest(data, engines = engines, horizons = horizons, ...)
  score_forecasts(bt, metrics = metrics)
}
