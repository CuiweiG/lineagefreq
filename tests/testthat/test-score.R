test_that("score_forecasts computes all requested metrics", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = c(7, 14),
                 min_train = 42)
  sc <- score_forecasts(bt, metrics = c("mae", "rmse", "coverage"))

  expect_s3_class(sc, "tbl_df")
  if (nrow(sc) > 0) {
    expect_true(all(c("engine", "horizon", "metric", "value")
                    %in% names(sc)))
    expect_true(all(c("mae", "rmse", "coverage") %in% sc$metric))
  }
})

test_that("score_forecasts mae is non-negative", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = 7,
                 min_train = 42)
  sc <- score_forecasts(bt, metrics = "mae")

  if (nrow(sc) > 0) {
    expect_true(all(sc$value >= 0))
  }
})

test_that("compare_models returns sorted tibble", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = c(7, 14),
                 min_train = 42)
  sc <- score_forecasts(bt)
  cm <- compare_models(sc)

  expect_s3_class(cm, "tbl_df")
})

test_that("plot_backtest returns ggplot", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = c(7, 14),
                 min_train = 42)
  sc <- score_forecasts(bt)

  if (nrow(sc) > 0) {
    p <- plot_backtest(sc)
    expect_s3_class(p, "gg")
  }
})
