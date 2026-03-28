# Helper: build multi-location lfq_data
make_hier_data <- function(seed1 = 1, seed2 = 2) {
  d1 <- simulate_dynamics(n_lineages = 3,
                          advantages = c("A" = 1.3, "B" = 0.8),
                          n_timepoints = 12, total_per_tp = 500, seed = seed1)
  d1_df <- as.data.frame(d1)
  d1_df$location <- "US"

  d2 <- simulate_dynamics(n_lineages = 3,
                          advantages = c("A" = 1.2, "B" = 0.9),
                          n_timepoints = 12, total_per_tp = 100, seed = seed2)
  d2_df <- as.data.frame(d2)
  d2_df$location <- "UK"

  combined <- rbind(d1_df, d2_df)
  lfq_data(combined, lineage = .lineage, date = .date,
           count = .count, location = location)
}

test_that("hier_mlr requires location column", {
  sim <- simulate_dynamics(n_lineages = 3,
                           advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
  expect_error(fit_model(sim, engine = "hier_mlr"), "location")
})

test_that("hier_mlr returns lfq_fit with correct class", {
  hdata <- make_hier_data()
  fit   <- fit_model(hdata, engine = "hier_mlr")

  expect_s3_class(fit, "lfq_fit")
  expect_s3_class(fit, "lfq_fit_hier_mlr")
  expect_equal(fit$engine, "hier_mlr")
})

test_that("hier_mlr has location-specific info", {
  hdata <- make_hier_data()
  fit   <- fit_model(hdata, engine = "hier_mlr")

  expect_true("locations" %in% names(fit))
  expect_equal(fit$n_locations, 2L)
  expect_true("location_fits" %in% names(fit))
  expect_equal(length(fit$location_fits), 2L)
})

test_that("hier_mlr produces shrinkage estimates", {
  hdata <- make_hier_data()
  fit   <- fit_model(hdata, engine = "hier_mlr")

  # Global growth rates should be numeric and named
  expect_true(is.numeric(fit$growth_rates))
  expect_true(all(c("A", "B", "ref") %in% names(fit$growth_rates)))

  # Shrinkage should exist if we have location fits
  if (!is.null(fit$shrunk_by_loc)) {
    expect_true(length(fit$shrunk_by_loc) > 0)
  }
})

test_that("hier_mlr S3 methods work", {
  hdata <- make_hier_data()
  fit   <- fit_model(hdata, engine = "hier_mlr")

  expect_output(print(fit))
  expect_true(is.numeric(coef(fit)))
  expect_true(is.matrix(vcov(fit)))

  ga <- growth_advantage(fit, type = "growth_rate")
  expect_s3_class(ga, "tbl_df")
  expect_equal(nrow(ga), 3L)
})
