test_that("MLR engine recovers known growth advantages", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("fast" = 1.4, "slow" = 0.7),
    n_timepoints = 20,
    total_per_tp = 2000,
    seed         = 42
  )
  # Explicitly set pivot = "ref" so the test is deterministic
  result <- .engine_mlr(sim, pivot = "ref")

  expect_true(is.list(result))
  expect_true("growth_rates" %in% names(result))

  # fast should have positive growth rate
  expect_gt(result$growth_rates["fast"], 0)
  # slow should have negative growth rate
  expect_lt(result$growth_rates["slow"], 0)
  # pivot (ref) should be exactly 0
  expect_equal(unname(result$growth_rates["ref"]), 0)
})

test_that("MLR engine auto-selects pivot", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.9),
    n_timepoints = 10, seed = 1
  )
  result <- .engine_mlr(sim)

  expect_true(result$pivot %in% attr(sim, "lineages"))
  expect_equal(unname(result$growth_rates[result$pivot]), 0)
})

test_that("MLR engine respects custom pivot", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.9),
    n_timepoints = 10, seed = 1
  )
  result <- .engine_mlr(sim, pivot = "A")

  expect_equal(result$pivot, "A")
  expect_equal(unname(result$growth_rates["A"]), 0)
})

test_that("MLR engine works with exactly 2 lineages", {
  sim <- simulate_dynamics(
    n_lineages   = 2,
    advantages   = c("rising" = 1.3),
    n_timepoints = 10,
    total_per_tp = 500,
    seed         = 5
  )
  result <- .engine_mlr(sim)

  expect_equal(length(result$lineages), 2)
  expect_equal(result$df, 2L)  # 1 alpha + 1 delta
})

test_that("MLR engine returns complete structure", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 12, seed = 7
  )
  result <- .engine_mlr(sim)

  required_fields <- c(
    "growth_rates", "intercepts", "pivot", "lineages",
    "fitted_values", "residuals", "vcov_matrix",
    "loglik", "aic", "bic", "nobs", "n_timepoints",
    "df", "convergence", "date_range", "time_scale"
  )
  for (f in required_fields) {
    expect_true(f %in% names(result),
                info = paste("Missing field:", f))
  }

  expect_true(is.matrix(result$vcov_matrix))
  expect_s3_class(result$fitted_values, "tbl_df")
  expect_s3_class(result$residuals, "tbl_df")
  expect_equal(result$convergence, 0L)
})

test_that("MLR engine window parameter reduces data", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 20, seed = 1
  )

  result_full   <- .engine_mlr(sim)
  result_window <- .engine_mlr(sim, window = 42)

  expect_lt(result_window$n_timepoints, result_full$n_timepoints)
})
