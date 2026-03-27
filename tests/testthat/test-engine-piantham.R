test_that("Piantham engine requires generation_time", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  expect_error(
    fit_model(sim, engine = "piantham"),
    "generation_time"
  )
})

test_that("Piantham engine returns lfq_fit with correct class", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  fit <- fit_model(sim, engine = "piantham", generation_time = 5)

  expect_s3_class(fit, "lfq_fit")
  expect_s3_class(fit, "lfq_fit_piantham")
  expect_equal(fit$engine, "piantham")
})

test_that("Piantham relative_Rt is consistent with growth_advantage", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.3, "B" = 0.8),
                           n_timepoints = 15, seed = 42)
  fit <- fit_model(sim, engine = "piantham", generation_time = 5)

  expect_true("relative_Rt" %in% names(fit))
  expect_equal(length(fit$relative_Rt), 3L)
  expect_equal(unname(fit$relative_Rt[fit$pivot]), 1)

  ga <- growth_advantage(fit, type = "relative_Rt", generation_time = 5)
  for (lin in fit$lineages) {
    expected_rt <- unname(fit$relative_Rt[lin])
    ga_rt       <- unname(ga$estimate[ga$lineage == lin])
    expect_equal(expected_rt, ga_rt, tolerance = 1e-6)
  }
})

test_that("Piantham works with all S3 methods", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  fit <- fit_model(sim, engine = "piantham", generation_time = 5)

  expect_output(print(fit))
  expect_true(is.numeric(coef(fit)))
  expect_s3_class(tidy.lfq_fit(fit), "tbl_df")
  expect_equal(nrow(glance.lfq_fit(fit)), 1L)
  expect_equal(glance.lfq_fit(fit)$engine, "piantham")
})
