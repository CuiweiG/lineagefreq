# Tests for joint conformal prediction on the simplex

test_that("ILR transform produces K-1 coordinates", {
  x <- c(0.5, 0.3, 0.2)
  ilr <- .ilr_transform(x)
  expect_length(ilr, 2)
  expect_true(all(is.finite(ilr)))
})

test_that("ILR transform is invertible", {
  x <- c(0.5, 0.3, 0.15, 0.05)
  ilr <- .ilr_transform(x)
  x_back <- .ilr_inverse(ilr)
  expect_equal(sum(x_back), 1, tolerance = 1e-10)
  expect_equal(x_back, x, tolerance = 1e-6)
})

test_that("Aitchison distance is zero for identical compositions", {
  x <- c(0.4, 0.35, 0.25)
  expect_equal(.aitchison_distance(x, x), 0, tolerance = 1e-10)
})

test_that("Aitchison distance is symmetric", {
  x <- c(0.5, 0.3, 0.2)
  y <- c(0.2, 0.5, 0.3)
  expect_equal(.aitchison_distance(x, y), .aitchison_distance(y, x),
               tolerance = 1e-10)
})

test_that("Aitchison distance is positive for different compositions", {
  x <- c(0.5, 0.3, 0.2)
  y <- c(0.2, 0.5, 0.3)
  expect_gt(.aitchison_distance(x, y), 0)
})

test_that("projected Aitchison ball respects simplex", {
  centre <- c(a = 0.5, b = 0.3, c = 0.2)
  bounds <- .project_aitchison_ball(centre, radius = 0.5,
                                    lineage_names = c("a", "b", "c"))
  # All bounds in [0, 1]
  expect_true(all(bounds$lower >= 0))
  expect_true(all(bounds$upper <= 1))
  # Lower <= centre <= upper
  expect_true(all(bounds$lower <= centre + 1e-6))
  expect_true(all(bounds$upper >= centre - 1e-6))
})

test_that("projected ball widens with larger radius", {
  centre <- c(a = 0.5, b = 0.3, c = 0.2)
  narrow <- .project_aitchison_ball(centre, radius = 0.1,
                                    lineage_names = c("a", "b", "c"))
  wide   <- .project_aitchison_ball(centre, radius = 1.0,
                                    lineage_names = c("a", "b", "c"))
  narrow_width <- narrow$upper - narrow$lower
  wide_width   <- wide$upper - wide$lower
  expect_true(all(wide_width >= narrow_width - 1e-6))
})

test_that("conformal_forecast_joint runs on simulated data", {
  skip_if_not_installed("lineagefreq")

  # Use built-in simulated data
  dat <- lineagefreq::sarscov2_us_2022
  lfq <- lfq_data(dat, lineage = variant, date = date, count = count)
  fit <- fit_model(lfq, engine = "mlr")

  result <- conformal_forecast_joint(fit, lfq, horizon = 14L,
                                      ci_level = 0.95,
                                      cal_fraction = 0.3,
                                      seed = 42)

  expect_s3_class(result, "lfq_conformal_joint")
  expect_gt(result$radius, 0)
  expect_gt(result$n_cal, 0)
  expect_true(nrow(result$marginal_intervals) > 0)

  # Joint intervals should have lower <= median <= upper
  mi <- result$marginal_intervals
  expect_true(all(mi$.lower_joint <= mi$.median + 1e-6))
  expect_true(all(mi$.upper_joint >= mi$.median - 1e-6))
})

test_that("calibrate_joint runs on backtest", {
  skip_if_not_installed("lineagefreq")

  dat <- lineagefreq::sarscov2_us_2022
  lfq <- lfq_data(dat, lineage = variant, date = date, count = count)
  bt  <- backtest(lfq, engines = "mlr", horizons = 14L, min_train = 42L)

  report <- calibrate_joint(bt)

  expect_s3_class(report, "joint_calibration_report")
  expect_gt(report$mean_energy_score, 0)
  expect_true(nrow(report$joint_coverage) > 0)
  expect_true(nrow(report$rank_histogram) > 0)
  expect_equal(sum(report$rank_histogram$count), report$n)
})

test_that("joint coverage is monotonically non-decreasing", {
  skip_if_not_installed("lineagefreq")

  dat <- lineagefreq::sarscov2_us_2022
  lfq <- lfq_data(dat, lineage = variant, date = date, count = count)
  bt  <- backtest(lfq, engines = "mlr", horizons = 14L, min_train = 42L)

  report <- calibrate_joint(bt)
  cov <- report$joint_coverage$observed
  # Coverage at higher nominal levels should be >= coverage at lower levels
  expect_true(all(diff(cov) >= -1e-10))
})
