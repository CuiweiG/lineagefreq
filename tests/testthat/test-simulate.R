test_that("simulate_dynamics returns valid lfq_data", {
  sim <- simulate_dynamics(seed = 42)

  expect_s3_class(sim, "lfq_data")
  expect_equal(length(attr(sim, "lineages")), 4)
  expect_equal(attr(sim, "n_timepoints"), 20L)
  expect_true(all(sim$.count >= 0))
})

test_that("named advantages become lineage names", {
  sim <- simulate_dynamics(
    n_lineages = 3,
    advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
    seed = 1
  )

  expect_true("JN.1" %in% attr(sim, "lineages"))
  expect_true("KP.3" %in% attr(sim, "lineages"))
  expect_true("ref"  %in% attr(sim, "lineages"))
})

test_that("seed gives reproducibility", {
  s1 <- simulate_dynamics(seed = 42)
  s2 <- simulate_dynamics(seed = 42)
  expect_equal(s1$.count, s2$.count)
})

test_that("wrong advantages length errors", {
  expect_error(
    simulate_dynamics(n_lineages = 3, advantages = c(1.2)),
    "length"
  )
})

test_that("overdispersion produces valid output", {
  sim <- simulate_dynamics(
    n_lineages     = 3,
    advantages     = c("A" = 1.2, "B" = 0.8),
    overdispersion = 0.05,
    seed           = 99
  )
  expect_s3_class(sim, "lfq_data")
  expect_true(all(sim$.count >= 0))
})

test_that("variable totals work correctly", {
  totals <- c(100L, 200L, 300L, 400L, 500L)
  sim <- simulate_dynamics(
    n_lineages   = 2,
    n_timepoints = 5,
    total_per_tp = totals,
    advantages   = c("A" = 1.1),
    seed         = 1
  )
  date_totals <- tapply(sim$.count, sim$.date, sum)
  expect_equal(as.integer(date_totals), totals)
})
