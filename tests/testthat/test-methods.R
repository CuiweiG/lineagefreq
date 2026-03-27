setup_fit <- function() {
  sim <- simulate_dynamics(
    n_lineages   = 3,
    advantages   = c("A" = 1.2, "B" = 0.8),
    n_timepoints = 12, seed = 7
  )
  fit_model(sim, engine = "mlr")
}

test_that("print.lfq_fit runs without error", {
  fit <- setup_fit()
  expect_output(print(fit), "Lineage frequency model")
})

test_that("summary.lfq_fit runs without error", {
  fit <- setup_fit()
  expect_output(summary(fit), "Model Summary")
})

test_that("coef returns named numeric", {
  fit <- setup_fit()

  gr <- coef(fit)
  expect_type(gr, "double")
  expect_equal(length(gr), 3L)   # 3 lineages
  expect_true(all(names(gr) %in% fit$lineages))

  all_coef <- coef(fit, type = "all")
  expect_equal(length(all_coef), 6L)  # 3 intercepts + 3 growth rates
})

test_that("vcov returns square matrix with correct names", {
  fit <- setup_fit()

  v <- vcov(fit)
  expect_true(is.matrix(v))
  expect_equal(nrow(v), ncol(v))
  expect_true(all(grepl("^(alpha|delta)_", colnames(v))))
})

test_that("logLik returns logLik class", {
  fit <- setup_fit()

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.numeric(ll))
  expect_equal(attr(ll, "df"), fit$df)
})

test_that("nobs returns integer", {
  fit <- setup_fit()
  expect_true(is.numeric(nobs(fit)))
  expect_gt(nobs(fit), 0)
})

test_that("tidy returns tibble with correct columns", {
  fit <- setup_fit()

  td <- tidy.lfq_fit(fit)
  expect_s3_class(td, "tbl_df")
  expect_true(all(c("lineage", "term", "estimate", "std.error",
                    "conf.low", "conf.high") %in% names(td)))
  # 2 non-pivot lineages x 2 terms = 4 rows
  expect_equal(nrow(td), 4L)
  expect_true(all(is.finite(td$estimate)))
})

test_that("tidy conf.int = FALSE omits intervals", {
  fit <- setup_fit()
  td  <- tidy.lfq_fit(fit, conf.int = FALSE)
  expect_true(all(is.na(td$conf.low)))
  expect_true(all(is.na(td$conf.high)))
})

test_that("glance returns single-row tibble", {
  fit <- setup_fit()

  gl <- glance.lfq_fit(fit)
  expect_s3_class(gl, "tbl_df")
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("engine", "n_lineages", "AIC", "BIC",
                    "logLik", "pivot") %in% names(gl)))
  expect_equal(gl$engine, "mlr")
})

test_that("augment returns residuals tibble", {
  fit <- setup_fit()

  aug <- augment.lfq_fit(fit)
  expect_s3_class(aug, "tbl_df")
  expect_true(all(c(".date", ".lineage", ".fitted_freq",
                    ".observed", ".pearson_resid") %in% names(aug)))
  expect_gt(nrow(aug), 0L)
})
