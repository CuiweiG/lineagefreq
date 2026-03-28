# v0.2.0 feature tests

test_that("pipe API works end-to-end", {
  data(sarscov2_us_2022)
  result <- sarscov2_us_2022 |>
    lfq_data(lineage = variant, date = date, count = count, total = total) |>
    lfq_fit("mlr") |>
    lfq_advantage(generation_time = 5)
  expect_s3_class(result, "tbl_df")
  expect_true("estimate" %in% names(result))
})

test_that("lfq_forecast pipe works", {
  sim <- simulate_dynamics(n_lineages = 2, advantages = c(A = 1.2), seed = 1)
  fc <- sim |> lfq_fit("mlr") |> lfq_forecast(14)
  expect_s3_class(fc, "lfq_forecast")
})

test_that("register_engine works", {
  my_fn <- function(data, pivot = NULL, ci_level = 0.95, ...) {
    .engine_mlr(data, pivot = pivot, ci_level = ci_level, ...)
  }
  register_engine("test_engine", my_fn, "Test wrapper")

  eng <- lfq_engines()
  expect_true("test_engine" %in% eng$engine)

  sim <- simulate_dynamics(n_lineages = 2, advantages = c(A = 1.1), seed = 1)
  fit <- fit_model(sim, engine = "test_engine")
  expect_s3_class(fit, "lfq_fit")

  unregister_engine("test_engine")
  eng2 <- lfq_engines()
  expect_false("test_engine" %in% eng2$engine)
})

test_that("lfq_summary works", {
  sim <- simulate_dynamics(n_lineages = 3, advantages = c(A = 1.2, B = 0.8), seed = 1)
  fit <- fit_model(sim)
  s <- lfq_summary(fit, generation_time = 5)
  expect_s3_class(s, "tbl_df")
  expect_true(all(c("lineage", "growth_rate", "relative_Rt") %in% names(s)))
})
