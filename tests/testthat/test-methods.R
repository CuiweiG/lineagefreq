make_fit <- function() {
  sim <- simulate_dynamics(
    n_lineages = 3, advantages = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 12, seed = 7
  )
  fit_model(sim, engine = "mlr")
}

test_that("print.lfq_fit runs without error", {
  fit <- make_fit()
  expect_no_error(print(fit))
  expect_invisible(print(fit))
})

test_that("summary.lfq_fit runs without error", {
  fit <- make_fit()
  expect_no_error(summary(fit))
})

test_that("coef.lfq_fit returns named numeric vector", {
  fit <- make_fit()
  cf  <- coef(fit)
  expect_true(is.numeric(cf))
  expect_true(!is.null(names(cf)))
  # Should have alpha and delta for each non-pivot
  n_non_pivot <- length(fit$lineages) - 1L
  expect_equal(length(cf), 2L * n_non_pivot)
})

test_that("vcov.lfq_fit returns a square matrix", {
  fit <- make_fit()
  vc  <- vcov(fit)
  expect_true(is.matrix(vc))
  expect_equal(nrow(vc), ncol(vc))
})

test_that("logLik.lfq_fit returns logLik object", {
  fit <- make_fit()
  ll  <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.numeric(as.numeric(ll)))
})

test_that("nobs.lfq_fit returns positive integer", {
  fit <- make_fit()
  n   <- nobs(fit)
  expect_true(is.integer(n))
  expect_gt(n, 0L)
})

test_that("tidy.lfq_fit returns correct structure", {
  fit <- make_fit()
  td  <- tidy.lfq_fit(fit)
  expect_s3_class(td, "tbl_df")
  expect_true(nrow(td) > 0)
  expect_true(all(c("lineage", "term", "estimate",
                    "std.error", "conf.low", "conf.high") %in% names(td)))
  # One alpha + one delta per non-pivot lineage
  n_non_pivot <- length(fit$lineages) - 1L
  expect_equal(nrow(td), 2L * n_non_pivot)
})

test_that("glance.lfq_fit returns one-row tibble", {
  fit <- make_fit()
  gl  <- glance.lfq_fit(fit)
  expect_s3_class(gl, "tbl_df")
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("engine", "n_lineages", "nobs",
                    "logLik", "AIC", "BIC") %in% names(gl)))
})

test_that("augment.lfq_fit returns residuals tibble", {
  fit <- make_fit()
  aug <- augment.lfq_fit(fit)
  expect_s3_class(aug, "tbl_df")
  expect_true(all(c(".date", ".lineage",
                    ".fitted_freq", ".observed") %in% names(aug)))
  expect_gt(nrow(aug), 0L)
})
