test_that("collapse_lineages merges rare lineages", {
  sim <- simulate_dynamics(
    n_lineages   = 5,
    advantages   = c(A = 1.3, B = 1.1, C = 0.5, D = 0.3),
    n_timepoints = 15, total_per_tp = 1000, seed = 1
  )
  collapsed <- collapse_lineages(sim, min_freq = 0.20)

  expect_s3_class(collapsed, "lfq_data")
  expect_true(
    length(attr(collapsed, "lineages")) <= length(attr(sim, "lineages"))
  )
  # At least some lineages should have been collapsed
  expect_true("Other" %in% attr(collapsed, "lineages"))
})

test_that("collapse_lineages with custom mapping", {
  sim <- simulate_dynamics(
    n_lineages   = 4,
    advantages   = c(A = 1.2, B = 0.9, C = 0.8),
    n_timepoints = 10, seed = 1
  )
  mapping   <- c("A" = "GroupX", "B" = "GroupX", "C" = "GroupY")
  collapsed <- collapse_lineages(sim, mapping = mapping)

  expect_true("GroupX" %in% attr(collapsed, "lineages"))
  expect_true("GroupY" %in% attr(collapsed, "lineages"))
})

test_that("filter_sparse removes low-count time points", {
  sim <- simulate_dynamics(
    n_lineages   = 2,
    advantages   = c("A" = 1.1),
    n_timepoints = 10,
    total_per_tp = c(5, 5, 500, 500, 500, 500, 500, 500, 500, 500),
    seed         = 1
  )
  filtered <- filter_sparse(sim, min_total = 100)

  expect_s3_class(filtered, "lfq_data")
  expect_true(all(filtered$.total >= 100))
})
