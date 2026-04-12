test_that("CRPS is non-negative and finite", {
  df <- tibble::tibble(
    predicted = 0.5,
    observed  = 0.5,
    lower     = 0.5 - 1.96 * 0.1,
    upper     = 0.5 + 1.96 * 0.1
  )
  crps <- .compute_crps(df)
  expect_true(is.finite(crps))
  expect_true(crps >= 0)

  # CRPS should increase when observation is far from prediction
  df2 <- tibble::tibble(
    predicted = 0.5,
    observed  = 0.9,
    lower     = 0.5 - 1.96 * 0.1,
    upper     = 0.5 + 1.96 * 0.1
  )
  crps2 <- .compute_crps(df2)
  expect_true(crps2 > crps)
})

test_that("log score is finite for reasonable forecasts", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = 7, min_train = 42)

  if (nrow(bt) > 0) {
    sc <- score_forecasts(bt, metrics = "log_score")
    expect_true(all(is.finite(sc$value)))
  }
})

test_that("DSS matches manual computation", {
  df <- tibble::tibble(
    predicted = 0.5,
    observed  = 0.6,
    lower     = 0.5 - 1.96 * 0.1,
    upper     = 0.5 + 1.96 * 0.1
  )
  dss <- .compute_dss(df)
  sigma <- 0.1
  expected <- log(sigma^2) + ((0.6 - 0.5) / sigma)^2
  expect_equal(dss, expected, tolerance = 0.01)
})

test_that("calibration score is zero for perfect calibration", {
  # If all observations land exactly at predicted, coverage should
  # be approximately correct at all levels
  n <- 200
  set.seed(42)
  sigma <- 0.1
  predicted <- rep(0.5, n)
  z <- stats::rnorm(n)
  observed <- predicted + z * sigma

  df <- tibble::tibble(
    predicted = predicted,
    observed  = observed,
    lower     = predicted - 1.96 * sigma,
    upper     = predicted + 1.96 * sigma
  )

  cal_score <- .compute_calibration_score(df)
  # With perfectly Gaussian data, calibration error should be small

  expect_true(cal_score < 0.01)
})

test_that("score_forecasts handles all new metrics", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = 7, min_train = 42)

  if (nrow(bt) > 0) {
    sc <- score_forecasts(bt, metrics = c("crps", "log_score", "dss",
                                          "calibration"))
    expect_true(all(c("crps", "log_score", "dss", "calibration")
                    %in% sc$metric))
    expect_true(all(is.finite(sc$value)))
  }
})

test_that("compare_models works with new scoring rules", {
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
    sc <- score_forecasts(bt, metrics = c("mae", "crps", "log_score"))
    cm <- compare_models(sc)
    expect_s3_class(cm, "tbl_df")
  }
})
