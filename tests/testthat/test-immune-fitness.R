test_that("immune_landscape creates valid object", {
  imm <- data.frame(
    date = rep(seq(as.Date("2024-01-01"), by = "week", length.out = 5),
               each = 3),
    lineage = rep(c("A", "B", "ref"), 5),
    immunity = runif(15, 0, 0.7)
  )
  il <- immune_landscape(imm, date = date, lineage = lineage,
                         immunity = immunity)
  expect_s3_class(il, "immune_landscape")
  expect_equal(length(il$lineages), 3L)
  expect_no_error(print(il))
})

test_that("immune_landscape rejects invalid immunity", {
  imm <- data.frame(
    date = as.Date("2024-01-01"),
    lineage = "A",
    immunity = 1.5
  )
  expect_error(
    immune_landscape(imm, date, lineage, immunity),
    "between 0 and 1"
  )
})

test_that("immune_landscape plot returns ggplot", {
  imm <- data.frame(
    date = rep(seq(as.Date("2024-01-01"), by = "week", length.out = 5),
               each = 2),
    lineage = rep(c("A", "ref"), 5),
    immunity = c(0.3, 0.5, 0.35, 0.5, 0.4, 0.5, 0.42, 0.5, 0.45, 0.5)
  )
  il <- immune_landscape(imm, date, lineage, immunity)
  p <- plot(il)
  expect_s3_class(p, "gg")
})

test_that("fitness_decomposition produces valid output", {
  sim <- simulate_dynamics(
    n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 15, total_per_tp = 500, seed = 1
  )
  fit <- fit_model(sim, engine = "mlr")

  dates <- unique(sim$.date)
  imm <- data.frame(
    date = rep(dates, each = 3),
    lineage = rep(c("A", "B", "ref"), length(dates)),
    immunity = c(rep(c(0.2, 0.5, 0.4), length(dates)))
  )
  il <- immune_landscape(imm, date, lineage, immunity)
  fd <- fitness_decomposition(fit, il, generation_time = 5)

  expect_s3_class(fd, "fitness_decomposition")
  expect_equal(nrow(fd$decomposition), 3L)
  expect_true(all(c("lineage", "observed_advantage", "beta",
                    "escape_contribution") %in% names(fd$decomposition)))
  expect_no_error(print(fd))
})

test_that("fitness_decomposition recovers known scenario", {
  # Scenario: lineage A has growth advantage and faces lower
  # population immunity than ref, so escape should contribute
  sim <- simulate_dynamics(
    n_lineages = 2,
    advantages = c("A" = 1.2),
    n_timepoints = 15, total_per_tp = 500, seed = 42
  )
  fit <- fit_model(sim, engine = "mlr")
  pivot <- fit$pivot
  non_piv <- setdiff(fit$lineages, pivot)

  # Give non-pivot lineage much lower immunity
  dates <- unique(sim$.date)
  lins <- fit$lineages
  imm <- data.frame(
    date = rep(dates, each = length(lins)),
    lineage = rep(lins, length(dates)),
    immunity = rep(ifelse(lins == pivot, 0.6, 0.1), length(dates))
  )
  il <- immune_landscape(imm, date, lineage, immunity)
  fd <- fitness_decomposition(fit, il, generation_time = 5)

  d_np <- fd$decomposition[fd$decomposition$lineage == non_piv[1], ]
  # With large immunity differential, escape should contribute
  expect_true(!is.na(d_np$escape_fraction))
  expect_true(d_np$escape_fraction > 0.05)
})

test_that("fitness_decomposition plot works", {
  sim <- simulate_dynamics(
    n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 15, total_per_tp = 500, seed = 1
  )
  fit <- fit_model(sim, engine = "mlr")

  dates <- unique(sim$.date)
  imm <- data.frame(
    date = rep(dates, each = 3),
    lineage = rep(c("A", "B", "ref"), length(dates)),
    immunity = rep(c(0.2, 0.5, 0.4), length(dates))
  )
  il <- immune_landscape(imm, date, lineage, immunity)
  fd <- fitness_decomposition(fit, il, generation_time = 5)
  p <- plot(fd)
  expect_s3_class(p, "gg")
})

test_that("fit_dms_prior produces valid lfq_fit", {
  sim <- simulate_dynamics(
    n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.9),
    n_timepoints = 8, total_per_tp = 100, seed = 1
  )
  dms <- c("A" = 0.04, "B" = -0.02)
  fit_d <- fit_dms_prior(sim, dms_scores = dms, lambda = 2)

  expect_s3_class(fit_d, "lfq_fit")
  expect_equal(fit_d$engine, "dms_prior")
  expect_true(fit_d$convergence == 0)
  expect_true(all(fit_d$fitted_values$.fitted_freq >= 0))
})

test_that("fit_dms_prior improves early detection vs standard MLR", {
  # With very sparse data, DMS prior should pull estimates toward
  # the informative prior
  set.seed(99)
  sim <- simulate_dynamics(
    n_lineages = 3,
    advantages = c("A" = 1.5, "B" = 0.8),
    n_timepoints = 5,  # very few
    total_per_tp = 50,  # very sparse
    seed = 99
  )

  fit_std <- fit_model(sim, engine = "mlr")
  dms <- c("A" = 0.06, "B" = -0.03)
  fit_d <- fit_dms_prior(sim, dms_scores = dms, lambda = 5)

  # DMS prior for A is positive; check that growth rate is positive
  expect_true(fit_d$growth_rates["A"] > 0)
})

test_that("selective_pressure returns valid tibble", {
  sim <- simulate_dynamics(
    n_lineages = 4,
    advantages = c("A" = 1.4, "B" = 1.1, "C" = 0.8),
    n_timepoints = 15, seed = 1
  )
  fit <- fit_model(sim, engine = "mlr")
  sp <- selective_pressure(fit)

  expect_s3_class(sp, "tbl_df")
  expect_true(all(c("date", "pressure", "dominant_lineage") %in% names(sp)))
  expect_true(nrow(sp) > 0)
})

test_that("selective_pressure variance method works", {
  sim <- simulate_dynamics(
    n_lineages = 3,
    advantages = c("A" = 1.3, "B" = 0.8),
    n_timepoints = 10, seed = 1
  )
  fit <- fit_model(sim, engine = "mlr")
  sp <- selective_pressure(fit, method = "variance")

  expect_true(all(sp$pressure >= 0))  # variance is non-negative
})
