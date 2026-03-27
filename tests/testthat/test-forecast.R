test_that("forecast returns lfq_forecast with correct structure", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  fit <- fit_model(sim, engine = "mlr")
  fc  <- forecast(fit, horizon = 14)

  expect_s3_class(fc, "lfq_forecast")
  expect_true("forecast" %in% fc$.type)
  expect_true("fitted"   %in% fc$.type)
  expect_true(all(c(".date", ".lineage", ".median", ".lower",
                    ".upper", ".type") %in% names(fc)))
})

test_that("forecast dates extend beyond fitted range", {
  sim <- simulate_dynamics(n_lineages = 2,
                           advantages = c("A" = 1.1),
                           n_timepoints = 10, seed = 1)
  fit <- fit_model(sim, engine = "mlr")
  fc  <- forecast(fit, horizon = 21)

  fc_dates <- unique(fc$.date[fc$.type == "forecast"])
  expect_true(all(fc_dates > max(fit$date_range)))
})

test_that("forecast frequencies are valid probabilities", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.3, "B" = 0.8),
                           n_timepoints = 15, seed = 42)
  fit <- fit_model(sim, engine = "mlr")
  fc  <- forecast(fit, horizon = 28, n_sim = 200)

  fc_part <- fc[fc$.type == "forecast", ]
  expect_true(all(fc_part$.median >= 0 & fc_part$.median <= 1))
  expect_true(all(fc_part$.lower  >= 0))
  expect_true(all(fc_part$.upper  <= 1 + 1e-9))
  expect_true(all(fc_part$.lower  <= fc_part$.median + 1e-9))
  expect_true(all(fc_part$.upper  >= fc_part$.median - 1e-9))
})

test_that("forecast CI gets wider with higher ci_level", {
  sim <- simulate_dynamics(n_lineages = 2,
                           advantages = c("A" = 1.2),
                           n_timepoints = 15, seed = 42)
  fit <- fit_model(sim, engine = "mlr")

  fc_95 <- forecast(fit, horizon = 14, ci_level = 0.95, n_sim = 500)
  fc_50 <- forecast(fit, horizon = 14, ci_level = 0.50, n_sim = 500)

  width_95 <- mean(fc_95$.upper[fc_95$.type == "forecast"] -
                     fc_95$.lower[fc_95$.type == "forecast"])
  width_50 <- mean(fc_50$.upper[fc_50$.type == "forecast"] -
                     fc_50$.lower[fc_50$.type == "forecast"])
  expect_gt(width_95, width_50)
})

test_that("forecast print works", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  fit <- fit_model(sim, engine = "mlr")
  fc  <- forecast(fit, horizon = 14)
  expect_no_error(print(fc))
})

test_that("autoplot.lfq_forecast returns ggplot", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  fit <- fit_model(sim, engine = "mlr")
  fc  <- forecast(fit, horizon = 14, n_sim = 200)

  p <- autoplot(fc)
  expect_s3_class(p, "gg")
})

test_that("autoplot.lfq_fit works for all types", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  fit <- fit_model(sim, engine = "mlr")

  expect_s3_class(autoplot(fit, type = "frequency"),  "gg")
  expect_s3_class(autoplot(fit, type = "trajectory"), "gg")
  expect_s3_class(autoplot(fit, type = "residuals"),  "gg")
  expect_s3_class(
    autoplot(fit, type = "advantage", generation_time = 5), "gg"
  )
})
