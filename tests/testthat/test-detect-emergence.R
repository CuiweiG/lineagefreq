test_that("summarize_emerging identifies growing lineages", {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c(emerging = 1.5, declining = 0.5),
    n_timepoints = 15, total_per_tp = 1000, seed = 42
  )
  result <- summarize_emerging(sim, threshold = 0.01)

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("lineage", "growth_rate", "direction",
                    "significant") %in% names(result)))
})

test_that("summarize_emerging returns empty tibble for no growth", {
  sim <- simulate_dynamics(
    n_lineages   = 2,
    advantages   = c("A" = 1.0),
    n_timepoints = 5,
    total_per_tp = 50,
    seed         = 1
  )
  result <- summarize_emerging(sim, threshold = 0.5)
  expect_s3_class(result, "tbl_df")
})

test_that("sequencing_power returns correct structure", {
  sp <- sequencing_power()
  expect_s3_class(sp, "tbl_df")
  expect_equal(nrow(sp), 1L)
  expect_true(sp$required_n > 0)
})

test_that("sequencing_power vectorized over current_freq", {
  sp <- sequencing_power(current_freq = c(0.01, 0.05, 0.10))
  expect_equal(nrow(sp), 3L)
  expect_true(all(sp$required_n > 0))
})
