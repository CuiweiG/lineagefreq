test_that("calibrate works on backtest results", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = c(7, 14),
                 min_train = 42)

  if (nrow(bt) > 0) {
    cal <- calibrate(bt)
    expect_s3_class(cal, "calibration_report")
    expect_true(length(cal$pit_values) > 0)
    expect_true(all(cal$pit_values >= 0 & cal$pit_values <= 1))
    expect_true(cal$ks_test$statistic >= 0)
    expect_true(cal$ks_test$p_value >= 0 & cal$ks_test$p_value <= 1)
    expect_equal(nrow(cal$reliability), 9L)
  }
})

test_that("calibrate print method works", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = 7, min_train = 42)

  if (nrow(bt) > 0) {
    cal <- calibrate(bt)
    expect_no_error(print(cal))
  }
})

test_that("calibrate plot returns ggplot", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = 7, min_train = 42)

  if (nrow(bt) > 0) {
    cal <- calibrate(bt)
    p <- plot(cal)
    expect_s3_class(p, "gg")
    p2 <- plot(cal, type = "pit")
    expect_s3_class(p2, "gg")
  }
})

test_that("recalibrate produces valid forecast", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  fit <- fit_model(sim, engine = "mlr")
  fc  <- forecast(fit, horizon = 14)
  bt  <- backtest(sim, engines = "mlr", horizons = 7, min_train = 42)

  if (nrow(bt) > 0) {
    fc_recal <- recalibrate(fc, bt, method = "isotonic")
    expect_s3_class(fc_recal, "lfq_forecast")
    fc_part <- fc_recal[fc_recal$.type == "forecast", ]
    expect_true(all(fc_part$.lower >= 0))
    expect_true(all(fc_part$.upper <= 1 + 1e-9))
  }
})

test_that("recalibrate platt method works", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  fit <- fit_model(sim, engine = "mlr")
  fc  <- forecast(fit, horizon = 14)
  bt  <- backtest(sim, engines = "mlr", horizons = 7, min_train = 42)

  if (nrow(bt) > 0) {
    fc_recal <- recalibrate(fc, bt, method = "platt")
    expect_s3_class(fc_recal, "lfq_forecast")
  }
})
