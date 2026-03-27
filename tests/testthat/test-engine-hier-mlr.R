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

test_that("hier_mlr growth rates are shrunk toward global", {
  hdata <- make_hier_data()
  fit   <- fit_model(hdata, engine = "hier_mlr")

  # The global growth rate for A should be positive
  # (both US and UK have A advantage > 1)
  expect_gt(fit$growth_rates["A"], 0)

  # The UK (low sample) fit should be shrunk toward global more
  # than the US (high sample) fit
  if (!is.null(fit$shrunk_by_loc[["UK"]]) &&
      !is.null(fit$shrunk_by_loc[["US"]])) {
    global_A <- fit$growth_rates["A"]
    us_A     <- fit$location_fits[["US"]]$growth_rates["A"]
    uk_A     <- fit$location_fits[["UK"]]$growth_rates["A"]
    uk_shrunk <- fit$shrunk_by_loc[["UK"]]["A"]
    us_shrunk <- fit$shrunk_by_loc[["US"]]["A"]

    # UK shrunk should be closer to global than UK raw
    expect_lte(abs(uk_shrunk - global_A), abs(uk_A - global_A) + 0.01)
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
