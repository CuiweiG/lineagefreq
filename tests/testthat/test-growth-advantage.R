test_that("growth_advantage returns correct structure", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  fit <- fit_model(sim, engine = "mlr")

  ga <- growth_advantage(fit, type = "growth_rate")
  expect_s3_class(ga, "tbl_df")
  expect_equal(nrow(ga), 3)
  expect_true(all(c("lineage", "estimate", "lower", "upper",
                    "type", "pivot") %in% names(ga)))
})

test_that("pivot lineage gets reference value", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  fit <- fit_model(sim, engine = "mlr")

  ga_gr     <- growth_advantage(fit, type = "growth_rate")
  pivot_row <- ga_gr[ga_gr$lineage == fit$pivot, ]
  expect_equal(unname(pivot_row$estimate), 0)

  ga_rt        <- growth_advantage(fit, type = "relative_Rt",
                                   generation_time = 5)
  pivot_row_rt <- ga_rt[ga_rt$lineage == fit$pivot, ]
  expect_equal(unname(pivot_row_rt$estimate), 1)
})

test_that("relative_Rt requires generation_time", {
  sim <- simulate_dynamics(n_lineages = 2,
                           advantages = c("A" = 1.2), seed = 1)
  fit <- fit_model(sim, engine = "mlr")

  expect_error(
    growth_advantage(fit, type = "relative_Rt"),
    "generation_time"
  )
})

test_that("all four types work without error", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.3, "B" = 0.8),
                           n_timepoints = 15, seed = 42)
  fit <- fit_model(sim, engine = "mlr")

  for (tp in c("growth_rate", "relative_Rt",
               "selection_coefficient", "doubling_time")) {
    ga <- growth_advantage(fit, type = tp, generation_time = 5)
    expect_s3_class(ga, "tbl_df")
    expect_equal(nrow(ga), 3)
    expect_true(all(is.finite(ga$estimate) | is.infinite(ga$estimate)))
  }
})

test_that("CI direction is correct for growth_rate", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8),
                           n_timepoints = 20, total_per_tp = 2000,
                           seed = 42)
  fit <- fit_model(sim, engine = "mlr")

  ga        <- growth_advantage(fit, type = "growth_rate")
  non_pivot <- ga[ga$lineage != fit$pivot, ]
  expect_true(all(non_pivot$lower <= non_pivot$estimate))
  expect_true(all(non_pivot$upper >= non_pivot$estimate))
})
