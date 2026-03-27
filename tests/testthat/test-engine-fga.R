test_that("FGA engine requires CmdStan", {
  skip_if(lfq_stan_available(),
          "Skipping: CmdStan IS available, testing error path not possible")

  sim <- simulate_dynamics(n_lineages = 2,
                           advantages = c("A" = 1.2), seed = 1)
  expect_error(fit_model(sim, engine = "fga"), "CmdStan")
})

test_that("FGA engine runs with CmdStan available", {
  skip_if_no_stan()

  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.3, "B" = 0.8),
    n_timepoints = 10,
    total_per_tp = 500,
    seed         = 42
  )

  fit <- fit_model(sim, engine = "fga",
                   chains = 2, iter_warmup = 200,
                   iter_sampling = 200)

  expect_s3_class(fit, "lfq_fit")
  expect_s3_class(fit, "lfq_fit_fga")
  expect_equal(fit$engine, "fga")
  expect_true(!is.null(fit$stan_fit))
  expect_true(!is.null(fit$rho_summary))
})

test_that("FGA growth_advantage works", {
  skip_if_no_stan()

  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.3, "B" = 0.8),
    n_timepoints = 10,
    total_per_tp = 500,
    seed         = 42
  )

  fit <- fit_model(sim, engine = "fga",
                   chains = 2, iter_warmup = 200,
                   iter_sampling = 200)

  ga <- growth_advantage(fit, type = "growth_rate")
  expect_s3_class(ga, "tbl_df")
  expect_equal(nrow(ga), 3L)
})

test_that("FGA S3 methods work", {
  skip_if_no_stan()

  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.3, "B" = 0.8),
    n_timepoints = 10,
    total_per_tp = 500,
    seed         = 42
  )

  fit <- fit_model(sim, engine = "fga",
                   chains = 2, iter_warmup = 200,
                   iter_sampling = 200)

  expect_no_error(print(fit))
  expect_true(is.numeric(coef(fit)))
  expect_s3_class(tidy.lfq_fit(fit), "tbl_df")
  expect_equal(nrow(glance.lfq_fit(fit)), 1L)
})
