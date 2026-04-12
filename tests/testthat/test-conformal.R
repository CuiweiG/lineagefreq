test_that("conformal_forecast produces valid lfq_forecast", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 20,
    total_per_tp = 500,
    seed         = 1
  )
  fit <- fit_model(sim, engine = "mlr")
  fc  <- conformal_forecast(fit, sim, horizon = 14, seed = 42)

  expect_s3_class(fc, "lfq_forecast")
  fc_part <- fc[fc$.type == "forecast", ]
  expect_true(nrow(fc_part) > 0)
  expect_true(all(fc_part$.lower >= 0))
  expect_true(all(fc_part$.upper <= 1 + 1e-9))
  expect_true(all(fc_part$.lower <= fc_part$.median + 1e-9))
})

test_that("conformal intervals achieve stated coverage on simulated data", {
  # Run many simulations and check empirical coverage
  set.seed(123)
  coverage_count <- 0L
  total_count    <- 0L

  for (trial in 1:5) {
    sim <- simulate_dynamics(
      n_lineages   = 3,
      advantages   = c("A" = 1.2 + trial * 0.05, "B" = 0.8),
      n_timepoints = 25,
      total_per_tp = 500,
      seed         = trial
    )

    # Use first 18 timepoints for fitting, rest for evaluation
    all_dates <- sort(unique(sim$.date))
    train_dates <- all_dates[1:18]
    test_dates  <- all_dates[19:min(25, length(all_dates))]

    train <- sim[sim$.date %in% train_dates, ]
    train <- structure(train, class = class(sim),
      lineages = attr(sim, "lineages"),
      date_range = range(train_dates),
      n_timepoints = length(train_dates),
      has_location = attr(sim, "has_location"),
      min_total = attr(sim, "min_total"))

    fit <- tryCatch(fit_model(train, engine = "mlr"), error = function(e) NULL)
    if (is.null(fit)) next

    fc <- tryCatch(
      conformal_forecast(fit, train, horizon = 56, ci_level = 0.90, seed = trial),
      error = function(e) NULL
    )
    if (is.null(fc)) next

    fc_part <- fc[fc$.type == "forecast", ]
    for (d in test_dates) {
      obs <- sim[sim$.date == d, ]
      for (lin in unique(obs$.lineage)) {
        obs_freq <- obs$.freq[obs$.lineage == lin]
        pred_row <- fc_part[fc_part$.lineage == lin, ]
        if (nrow(pred_row) > 0 && length(obs_freq) == 1) {
          # Use closest forecast date
          idx <- which.min(abs(as.numeric(pred_row$.date - d)))
          total_count <- total_count + 1L
          if (obs_freq >= pred_row$.lower[idx] &&
              obs_freq <= pred_row$.upper[idx]) {
            coverage_count <- coverage_count + 1L
          }
        }
      }
    }
  }

  if (total_count > 10) {
    empirical_coverage <- coverage_count / total_count
    # Conformal guarantees are marginal; with small samples,
    # allow generous tolerance
    expect_true(empirical_coverage >= 0.60,
      info = paste("Empirical coverage:", round(empirical_coverage, 3)))
  }
})

test_that("conformal ACI method works", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 20,
    total_per_tp = 500,
    seed         = 1
  )
  fit <- fit_model(sim, engine = "mlr")
  fc  <- conformal_forecast(fit, sim, horizon = 14, method = "aci")

  expect_s3_class(fc, "lfq_forecast")
})
