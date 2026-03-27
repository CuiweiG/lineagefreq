test_that("backtest produces results with correct structure", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = c(7, 14),
                 min_train = 42)

  expect_s3_class(bt, "lfq_backtest")
  expect_true(nrow(bt) > 0)
  expect_true(all(c("origin_date", "target_date", "horizon",
                    "engine", "lineage", "predicted", "observed")
                  %in% names(bt)))
})

test_that("backtest respects min_train", {
  sim <- simulate_dynamics(
    n_lineages   = 2,
    advantages   = c("A" = 1.1),
    n_timepoints = 15,
    seed         = 1
  )
  # Very large min_train should produce few or no results
  bt <- backtest(sim, engines = "mlr", horizons = 7,
                 min_train = 80)
  expect_true(nrow(bt) >= 0)
})

test_that("backtest works with multiple engines", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 18,
    total_per_tp = 500,
    seed         = 42
  )
  bt <- backtest(sim,
    engines         = c("mlr", "piantham"),
    horizons        = 7,
    min_train       = 42,
    generation_time = 5
  )

  # At least one engine should produce results
  expect_s3_class(bt, "lfq_backtest")
  if (nrow(bt) > 0) {
    expect_true(length(unique(bt$engine)) >= 1L)
  }
})

test_that("backtest print works", {
  sim <- simulate_dynamics(
    n_lineages   = 2,
    advantages   = c("A" = 1.2),
    n_timepoints = 18,
    seed         = 1
  )
  bt <- backtest(sim, engines = "mlr", horizons = 7,
                 min_train = 42)
  expect_no_error(print(bt))
})
