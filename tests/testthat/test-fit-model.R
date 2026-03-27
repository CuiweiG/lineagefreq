test_that("fit_model returns lfq_fit with correct class", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  fit <- fit_model(sim, engine = "mlr")

  expect_s3_class(fit, "lfq_fit")
  expect_s3_class(fit, "lfq_fit_mlr")
  expect_equal(fit$engine, "mlr")
  expect_true(!is.null(fit$call))
})

test_that("fit_model rejects non-lfq_data input", {
  d <- data.frame(x = 1)
  expect_error(fit_model(d, engine = "mlr"), "lfq_data")
})

test_that("fit_model hier_mlr placeholder gives informative error", {
  sim <- simulate_dynamics(n_lineages = 2,
                           advantages = c("A" = 1.1), seed = 1)
  expect_error(fit_model(sim, engine = "hier_mlr"), "not yet")
})

test_that("fit_model piantham requires generation_time", {
  sim <- simulate_dynamics(n_lineages = 2,
                           advantages = c("A" = 1.1), seed = 1)
  expect_error(fit_model(sim, engine = "piantham"), "generation_time")
})
