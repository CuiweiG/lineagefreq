# Tests for pseudo-prospective evaluation

test_that("evaluate_prospective runs on simulated data", {
  skip_if_not_installed("lineagefreq")

  dat <- lineagefreq::sarscov2_us_2022
  lfq <- lfq_data(dat, lineage = variant, date = date, count = count)

  result <- evaluate_prospective(lfq, engine = "mlr",
                                  horizons = 14L,
                                  min_train = 42L, min_cal = 3L,
                                  ci_level = 0.95, gamma = 0.05)

  expect_s3_class(result, "lfq_prospective")
  expect_gt(result$n_origins, 0)
  expect_gt(nrow(result$results), 0)
  expect_true(all(c("lower_param", "lower_static", "lower_aci") %in%
                    names(result$results)))
})

test_that("prospective results have valid structure", {
  skip_if_not_installed("lineagefreq")

  dat <- lineagefreq::sarscov2_us_2022
  lfq <- lfq_data(dat, lineage = variant, date = date, count = count)

  result <- evaluate_prospective(lfq, engine = "mlr",
                                  horizons = 14L,
                                  min_train = 42L, min_cal = 3L)

  res <- result$results
  # Predictions in [0, 1]
  expect_true(all(res$predicted >= 0 & res$predicted <= 1, na.rm = TRUE))
  expect_true(all(res$observed >= 0 & res$observed <= 1, na.rm = TRUE))

  # Static intervals widen or stay stable over time (calibration set grows)
  valid <- res[!is.na(res$radius_static), ]
  if (nrow(valid) > 2) {
    # n_cal should be non-decreasing
    expect_true(all(diff(valid$n_cal) >= 0))
  }
})

test_that("summary contains three methods", {
  skip_if_not_installed("lineagefreq")

  dat <- lineagefreq::sarscov2_us_2022
  lfq <- lfq_data(dat, lineage = variant, date = date, count = count)

  result <- evaluate_prospective(lfq, engine = "mlr",
                                  horizons = 14L,
                                  min_train = 42L, min_cal = 3L)

  expect_equal(nrow(result$summary), 3)
  expect_equal(sort(result$summary$method),
               c("ACI", "Parametric", "Static conformal"))
  # Coverage values in [0, 1]
  expect_true(all(result$summary$coverage >= 0 &
                    result$summary$coverage <= 1))
})

test_that("ACI alpha adapts over time", {
  skip_if_not_installed("lineagefreq")

  dat <- lineagefreq::sarscov2_us_2022
  lfq <- lfq_data(dat, lineage = variant, date = date, count = count)

  result <- evaluate_prospective(lfq, engine = "mlr",
                                  horizons = 14L,
                                  min_train = 42L, min_cal = 3L,
                                  gamma = 0.1)

  valid <- result$results[!is.na(result$results$alpha_aci), ]
  if (nrow(valid) > 5) {
    # Alpha should vary (not all identical)
    expect_gt(stats::sd(valid$alpha_aci), 0)
  }
})

test_that("plot methods produce ggplot objects", {
  skip_if_not_installed("lineagefreq")

  dat <- lineagefreq::sarscov2_us_2022
  lfq <- lfq_data(dat, lineage = variant, date = date, count = count)

  result <- evaluate_prospective(lfq, engine = "mlr",
                                  horizons = 14L,
                                  min_train = 42L, min_cal = 3L)

  p1 <- plot(result, type = "coverage")
  expect_s3_class(p1, "gg")

  p2 <- plot(result, type = "radius")
  expect_s3_class(p2, "gg")

  p3 <- plot(result, type = "comparison")
  expect_s3_class(p3, "gg")
})
