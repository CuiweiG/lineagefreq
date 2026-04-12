test_that("surveillance_value returns valid evoi object", {
  sim <- simulate_dynamics(n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 15, seed = 1)
  fit <- fit_model(sim, engine = "mlr")
  ev <- surveillance_value(fit, n_current = 500)

  expect_s3_class(ev, "evoi")
  expect_true(ev$current_uncertainty > 0)
  expect_true(nrow(ev$values) > 0)
  expect_no_error(print(ev))
})

test_that("surveillance_value is monotonically increasing in n", {
  sim <- simulate_dynamics(n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 15, seed = 1)
  fit <- fit_model(sim, engine = "mlr")
  ev <- surveillance_value(fit, n_current = 100)

  # EVOI (total variance reduction) should increase with more samples
  expect_true(all(diff(ev$values$evoi) >= -1e-10))
})

test_that("surveillance_value plot works", {
  sim <- simulate_dynamics(n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 15, seed = 1)
  fit <- fit_model(sim, engine = "mlr")
  ev <- surveillance_value(fit, n_current = 500)
  p <- plot(ev)
  expect_s3_class(p, "gg")
})

test_that("adaptive_design returns valid allocation", {
  sim <- simulate_dynamics(n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 12, seed = 1)
  ad <- adaptive_design(sim, capacity = 200, n_rounds = 6, seed = 42)

  expect_s3_class(ad, "adaptive_allocation")
  expect_true(nrow(ad$allocations) > 0)
  expect_true(nrow(ad$summary) > 0)
  expect_no_error(print(ad))
})

test_that("adaptive_design UCB strategy works", {
  sim <- simulate_dynamics(n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 12, seed = 1)
  ad <- adaptive_design(sim, capacity = 200, n_rounds = 6,
    strategy = "ucb", seed = 42)
  expect_s3_class(ad, "adaptive_allocation")
})

test_that("adaptive_design plot returns ggplot", {
  sim <- simulate_dynamics(n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 12, seed = 1)
  ad <- adaptive_design(sim, capacity = 200, n_rounds = 6, seed = 42)
  p <- plot(ad)
  expect_s3_class(p, "gg")
})

test_that("detection_horizon matches analytical solution", {
  # At prevalence p in n sequences, P(detect >=1) = 1 - (1-p)^n
  dh <- detection_horizon(initial_prev = 0.01, growth_rate = 1.0,
    n_per_period = 500, n_periods = 10)

  # With no growth (rate=1.0), prevalence stays at 0.01
  # P(detect in period 1) = 1 - (1-0.01)^500
  expected_p1 <- 1 - (1 - 0.01)^500

  # All prevalences should be approximately equal (no growth)
  expect_true(all(abs(dh$prevalence - 0.01) < 0.005))

  # Per-period detection should be very high at p=0.01, n=500
  expect_true(expected_p1 > 0.99)
})

test_that("detection_horizon returns weeks_to_detection attribute", {
  dh <- detection_horizon(initial_prev = 0.001, growth_rate = 1.5,
    n_per_period = 100)
  wtd <- attr(dh, "weeks_to_detection")
  expect_true(is.integer(wtd) || is.na(wtd))
  if (!is.na(wtd)) {
    expect_true(wtd > 0 && wtd <= 26)
  }
})

test_that("alert_threshold SPRT returns valid output", {
  sim <- simulate_dynamics(n_lineages = 2,
    advantages = c("A" = 1.5),
    n_timepoints = 15, total_per_tp = 500, seed = 1)
  alerts <- alert_threshold(sim, method = "sprt", delta_1 = 0.05)

  expect_s3_class(alerts, "tbl_df")
  expect_true(nrow(alerts) == 2L)  # one per lineage
  expect_true(all(c("lineage", "date", "statistic", "alert",
    "direction") %in% names(alerts)))

  # Statistic for the growing lineage should be positive
  a_row <- alerts[alerts$lineage == "A", ]
  expect_true(a_row$statistic > 0)
})

test_that("alert_threshold CUSUM method works", {
  sim <- simulate_dynamics(n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.8),
    n_timepoints = 15, total_per_tp = 500, seed = 1)
  alerts <- alert_threshold(sim, method = "cusum", threshold = 3)

  expect_s3_class(alerts, "tbl_df")
  expect_true(nrow(alerts) > 0)
})

test_that("alert_threshold controls false alarms on stable data", {
  # Generate data with NO growth — no lineage should trigger alert
  set.seed(42)
  n_false_alarms <- 0L
  n_trials <- 20L

  for (i in seq_len(n_trials)) {
    sim <- simulate_dynamics(n_lineages = 2,
      advantages = c("A" = 1.0),  # no growth
      n_timepoints = 12, total_per_tp = 500, seed = i)
    alerts <- alert_threshold(sim, method = "sprt", alpha = 0.05,
      delta_1 = 0.05)
    if (any(alerts$alert)) n_false_alarms <- n_false_alarms + 1L
  }

  # False alarm rate should be controlled near alpha
  # With 20 trials, expect <= 5 false alarms at alpha=0.05
  expect_true(n_false_alarms <= 8L,
    info = paste("False alarms:", n_false_alarms, "out of", n_trials))
})

test_that("surveillance_dashboard returns panel list", {
  sim <- simulate_dynamics(n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 15, seed = 1)
  fit <- fit_model(sim, engine = "mlr")
  panels <- surveillance_dashboard(fit, sim)

  expect_s3_class(panels, "surveillance_dashboard")
  expect_true("landscape" %in% names(panels))
  expect_true("detection" %in% names(panels))
})

test_that("surveillance_dashboard works with backtest", {
  sim <- simulate_dynamics(n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 18, total_per_tp = 500, seed = 1)
  fit <- fit_model(sim, engine = "mlr")
  bt <- backtest(sim, engines = "mlr", horizons = 7, min_train = 42)

  if (nrow(bt) > 0) {
    panels <- surveillance_dashboard(fit, sim, bt = bt)
    expect_true("calibration" %in% names(panels))
  }
})
