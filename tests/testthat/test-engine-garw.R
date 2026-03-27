test_that("GARW engine requires CmdStan", {
  skip_if(lfq_stan_available(),
          "Skipping: CmdStan IS available")

  sim <- simulate_dynamics(n_lineages = 2,
                           advantages = c("A" = 1.2), seed = 1)
  expect_error(fit_model(sim, engine = "garw"), "CmdStan")
})

test_that("GARW engine runs with CmdStan available", {
  skip_if_no_stan()

  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.3, "B" = 0.8),
    n_timepoints = 10,
    total_per_tp = 500,
    seed         = 42
  )

  fit <- fit_model(sim, engine = "garw",
                   chains = 2, iter_warmup = 200,
                   iter_sampling = 200)

  expect_s3_class(fit, "lfq_fit")
  expect_s3_class(fit, "lfq_fit_garw")
  expect_equal(fit$engine, "garw")
  expect_true(fit$is_time_varying)
  expect_true(!is.null(fit$rho_trajectory))
  expect_true(!is.null(fit$sigma_rw))
})

test_that("GARW rho_trajectory has correct structure", {
  skip_if_no_stan()

  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.3, "B" = 0.8),
    n_timepoints = 10,
    total_per_tp = 500,
    seed         = 42
  )

  fit <- fit_model(sim, engine = "garw",
                   chains = 2, iter_warmup = 200,
                   iter_sampling = 200)

  rt <- fit$rho_trajectory
  expect_s3_class(rt, "tbl_df")
  expect_true(all(c(".date", ".lineage", "rho_median",
                    "rho_lower", "rho_upper") %in% names(rt)))
  expect_equal(nrow(rt), 10L * 3L)
})

test_that("GARW S3 methods work", {
  skip_if_no_stan()

  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.3, "B" = 0.8),
    n_timepoints = 10,
    total_per_tp = 500,
    seed         = 42
  )

  fit <- fit_model(sim, engine = "garw",
                   chains = 2, iter_warmup = 200,
                   iter_sampling = 200)

  expect_no_error(print(fit))
  expect_true(is.numeric(coef(fit)))
  expect_s3_class(tidy.lfq_fit(fit), "tbl_df")
  gl <- glance.lfq_fit(fit)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$engine, "garw")
})
