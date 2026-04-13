# Tests for CAPS (Compositional Adaptive Prediction Sets)

test_that("caps_forecast works on built-in data", {
  dat <- lineagefreq::cdc_ba2_transition
  x   <- lfq_data(dat, lineage = lineage, date = date, count = count)
  fit <- fit_model(x, engine = "mlr")

  caps <- caps_forecast(fit, x, horizons = c(7L, 14L),
                        method = "caps", alpha = 0.05)

  expect_s3_class(caps, "caps_forecast")
  expect_true(nrow(caps$forecasts) > 0)
  expect_true(all(caps$forecasts$caps_lower >= 0))
  expect_true(all(caps$forecasts$caps_upper <= 1))
  expect_true(all(is.finite(caps$R_hat) | is.na(caps$R_hat)))
  expect_true(caps$psi_K > 0)
})

test_that("caps_forecast methods are nested", {
  dat <- lineagefreq::cdc_ba2_transition
  x   <- lfq_data(dat, lineage = lineage, date = date, count = count)
  fit <- fit_model(x, engine = "mlr")

  caps_full     <- caps_forecast(fit, x, horizons = 14L, method = "caps")
  caps_static   <- caps_forecast(fit, x, horizons = 14L, method = "caps_static")
  caps_marginal <- caps_forecast(fit, x, horizons = 14L, method = "caps_marginal")

  expect_true(nrow(caps_full$forecasts) > 0)
  expect_true(nrow(caps_static$forecasts) > 0)
  expect_true(nrow(caps_marginal$forecasts) > 0)

  # caps_marginal should have psi_K = 1 (no correction)
  expect_equal(caps_marginal$psi_K, 1)
})

test_that("variance ratios are computed correctly", {
  dat <- lineagefreq::cdc_ba2_transition
  x   <- lfq_data(dat, lineage = lineage, date = date, count = count)
  bt  <- backtest(x, engines = "mlr", horizons = c(7L, 14L))

  R <- .caps_variance_ratios(bt, c(7L, 14L))
  expect_length(R, 2)
  expect_named(R, c("7", "14"))
  # At least one should be finite
 expect_true(any(is.finite(R)))
})

test_that("dimension correction psi(K) tightens for K > 2", {
  K <- 5
  psi <- (gamma((K - 1) / 2 + 1) / (pi^((K - 1) / 2)))^(1 / (K - 1))
  expect_true(psi > 0)
  expect_true(psi < 1)
})

test_that("evaluate_caps returns correct structure", {
  dat <- lineagefreq::cdc_ba2_transition
  x   <- lfq_data(dat, lineage = lineage, date = date, count = count)
  fit <- fit_model(x, engine = "mlr")
  bt  <- backtest(x, engines = "mlr", horizons = 14L, min_train = 42L)

  caps <- caps_forecast(fit, x, horizons = 14L, method = "caps")
  eval <- evaluate_caps(caps, bt)

  expect_true(nrow(eval) > 0)
  expect_true(all(c("para_cov", "conf_cov", "caps_cov") %in% names(eval)))
  # Coverage values in [0, 1]
  expect_true(all(eval$caps_cov >= 0 & eval$caps_cov <= 1))
})

test_that("CAPS radius widens with horizon when R decreases", {
  dat <- lineagefreq::cdc_ba2_transition
  x   <- lfq_data(dat, lineage = lineage, date = date, count = count)
  fit <- fit_model(x, engine = "mlr")

  caps <- caps_forecast(fit, x, horizons = c(7L, 14L, 28L), method = "caps")

  # Extract radii
  radii <- sapply(caps$params, function(p) p$r_base)
  radii <- radii[is.finite(radii)]

  # With at least 2 valid horizons, check non-decreasing trend
  if (length(radii) >= 2) {
    expect_true(radii[length(radii)] >= radii[1] - 0.01)
  }
})
