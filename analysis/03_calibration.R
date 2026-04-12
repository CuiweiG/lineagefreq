###############################################################################
# 03_calibration.R — Calibration diagnostics (Conclusions 2 + 3)
# lineagefreq validation analysis
###############################################################################

cat("[03_calibration] Starting...\n")
source("analysis/00_setup.R")

set.seed(20240413)

benchmark <- readRDS("analysis/results/benchmark_multicountry.rds")
cat(sprintf("  Loaded benchmark results: %d successful configurations\n",
            length(benchmark$backtest_results)))

###############################################################################
# Section A: PIT diagnostics across all countries
###############################################################################

cat("  [A] PIT diagnostics...\n")

# For each dataset-engine backtest, compute PIT values and KS test
# KEY HYPOTHESIS: If U-shape appears in ALL 5 countries + US,
# this proves miscalibration is inherent to MLR, not a data-specific bug.

pit_results <- list()

for (nm in names(benchmark$backtest_results)) {
  bt <- benchmark$backtest_results[[nm]]
  # backtest() returns an lfq_backtest tibble directly (no $forecasts)
  bt_tibble <- bt$result

  cat(sprintf("    Calibrating %s/%s...\n", bt$dataset, bt$engine))

  # Check if backtest has non-NA lower/upper values needed for calibration
  has_intervals <- all(c("lower", "upper") %in% names(bt_tibble)) &&
    any(!is.na(bt_tibble$lower)) && any(!is.na(bt_tibble$upper))

  if (!has_intervals) {
    cat(sprintf("      %s engine does not produce prediction intervals; calibration skipped\n",
                bt$engine))
    next
  }

  tryCatch({
    cal <- calibrate(bt_tibble)

    pit_results[[nm]] <- list(
      dataset    = bt$dataset,
      engine     = bt$engine,
      pit_values = cal$pit_values,
      ks_D       = cal$ks_test$statistic,
      ks_p       = cal$ks_test$p_value,
      n_pit      = cal$n,
      calibration_obj = cal
    )

    cat(sprintf("      KS D = %.3f, p = %.3g (n = %d)\n",
                cal$ks_test$statistic, cal$ks_test$p_value, cal$n))
  },
  error = function(e) {
    cat(sprintf("    WARNING: calibrate() failed for %s: %s\n", nm, e$message))
  })
}

cat(sprintf("    PIT diagnostics completed for %d configurations\n",
            length(pit_results)))

###############################################################################
# Section B: Reliability diagram data
###############################################################################

cat("  [B] Reliability diagram...\n")

# calibrate() already computes reliability (nominal vs observed coverage).
# Extract from calibration reports, plus compute for additional nominal levels.

nominal_levels <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

reliability_data <- list()

for (nm in names(pit_results)) {
  pr <- pit_results[[nm]]
  cal <- pr$calibration_obj

  tryCatch({
    # calibrate() returns reliability at 0.1-0.9 by 0.1
    # We also need 0.95; compute from PIT values
    pit <- pr$pit_values

    obs_coverage <- sapply(nominal_levels, function(level) {
      alpha <- (1 - level) / 2
      mean(pit >= alpha & pit <= 1 - alpha, na.rm = TRUE)
    })

    reliability_data[[nm]] <- tibble(
      dataset           = pr$dataset,
      engine            = pr$engine,
      nominal_coverage  = nominal_levels,
      observed_coverage = obs_coverage
    )
  },
  error = function(e) {
    cat(sprintf("    WARNING: Reliability data failed for %s: %s\n", nm, e$message))
  })
}

reliability_df <- bind_rows(reliability_data)
cat(sprintf("    Reliability data for %d configurations\n", length(reliability_data)))

###############################################################################
# Section C: Conformal vs Parametric vs Recalibrated comparison
###############################################################################

cat("  [C] Conformal vs Parametric vs Recalibrated...\n")

# ACI theoretical guarantee:
# Gibbs & Candès (2021, "Adaptive Conformal Inference Under Distribution Shift",
# NeurIPS) prove that adaptive conformal inference achieves asymptotic marginal
# coverage guarantees. For non-exchangeable data (as in time series), the coverage
# guarantee is asymptotic and depends on the step-size parameter γ.
# In finite samples with strong distribution shift (variant transitions), coverage
# may deviate from nominal. This is a known limitation.

# conformal_forecast() and recalibrate() require lfq_forecast/lfq_fit objects,
# not backtest tibbles. Instead, we compute conformal intervals from backtest
# residuals directly using split conformal prediction:
#   1. Split backtest origins into calibration (first 60%) and test (last 40%)
#   2. Compute absolute residuals on calibration set
#   3. Use the (1-α)(1+1/n_cal) quantile of absolute residuals as conformal radius
#   4. Apply to test set predictions: [predicted ± radius]
# This is valid split conformal prediction (Lei et al. 2018, JRSS-B).

# Filter to MLR engine (piantham lacks intervals)
mlr_backtests <- Filter(function(x) x$engine == "mlr", benchmark$backtest_results)

calibration_comparison <- list()

for (nm in names(mlr_backtests)) {
  bt <- mlr_backtests[[nm]]
  bt_tibble <- bt$result

  cat(sprintf("    Comparing methods for %s...\n", bt$dataset))

  tryCatch({
    # Split by origin_date: first 60% calibration, last 40% test
    origin_dates <- sort(unique(bt_tibble$origin_date))
    n_origins    <- length(origin_dates)
    n_cal        <- floor(n_origins * 0.6)

    if (n_cal < 3 || n_origins - n_cal < 3) {
      cat("      WARNING: Too few origins for split conformal; skipping\n")
      next
    }

    cal_origins  <- origin_dates[1:n_cal]
    test_origins <- origin_dates[(n_cal + 1):n_origins]

    cal_set  <- bt_tibble |> filter(origin_date %in% cal_origins)
    test_set <- bt_tibble |> filter(origin_date %in% test_origins)

    # ── 1. Parametric (existing lower/upper from backtest) ────────────────
    parametric_test <- test_set |>
      mutate(method = "Parametric")

    para_cov95  <- mean(test_set$observed >= test_set$lower &
                         test_set$observed <= test_set$upper, na.rm = TRUE)
    para_width  <- mean(test_set$upper - test_set$lower, na.rm = TRUE)

    # Winkler score (interval score) at 95%
    # Winkler (1972): S = width + (2/α)(lower - obs)[obs<lower]
    #                            + (2/α)(obs - upper)[obs>upper]
    alpha <- 0.05
    para_winkler <- mean(
      (test_set$upper - test_set$lower) +
        (2 / alpha) * pmax(0, test_set$lower - test_set$observed) +
        (2 / alpha) * pmax(0, test_set$observed - test_set$upper),
      na.rm = TRUE
    )

    # ── 2. Conformal (split conformal from residuals) ─────────────────────
    cal_residuals <- abs(cal_set$predicted - cal_set$observed)
    # Conformal quantile: (1 - α)(1 + 1/n_cal) quantile
    # (Lei et al. 2018; guarantees marginal coverage for exchangeable data)
    conf_level   <- (1 - alpha) * (1 + 1 / length(cal_residuals))
    conf_level   <- min(conf_level, 1)  # clip to [0, 1]
    conf_radius  <- quantile(cal_residuals, probs = conf_level, na.rm = TRUE)

    conformal_test <- test_set |>
      mutate(
        method       = "Conformal",
        conf_lower   = predicted - conf_radius,
        conf_upper   = predicted + conf_radius
      )

    conf_cov95   <- mean(test_set$observed >= conformal_test$conf_lower &
                          test_set$observed <= conformal_test$conf_upper,
                         na.rm = TRUE)
    conf_width   <- 2 * conf_radius
    conf_winkler <- mean(
      conf_width +
        (2 / alpha) * pmax(0, conformal_test$conf_lower - test_set$observed) +
        (2 / alpha) * pmax(0, test_set$observed - conformal_test$conf_upper),
      na.rm = TRUE
    )

    # ── 3. Recalibrated (scale parametric intervals by calibration ratio) ─
    # Estimate the needed scaling factor from calibration set
    cal_covered <- cal_set$observed >= cal_set$lower &
                   cal_set$observed <= cal_set$upper
    cal_coverage_obs <- mean(cal_covered, na.rm = TRUE)

    # If parametric intervals undercover, widen by a factor
    # Find scaling factor s such that coverage of [pred ± s*halfwidth] ≈ 95%
    cal_halfwidth <- (cal_set$upper - cal_set$lower) / 2
    cal_z_scores  <- abs(cal_set$predicted - cal_set$observed) /
                     pmax(cal_halfwidth, 1e-10)

    # The 95th percentile of |z| gives the needed scaling factor
    recal_factor <- quantile(cal_z_scores, probs = 0.95, na.rm = TRUE)
    recal_factor <- max(recal_factor, 1)  # never shrink intervals

    test_halfwidth <- (test_set$upper - test_set$lower) / 2
    recalibrated_test <- test_set |>
      mutate(
        method       = "Recalibrated",
        recal_lower  = predicted - recal_factor * test_halfwidth,
        recal_upper  = predicted + recal_factor * test_halfwidth
      )

    recal_cov95  <- mean(test_set$observed >= recalibrated_test$recal_lower &
                          test_set$observed <= recalibrated_test$recal_upper,
                         na.rm = TRUE)
    recal_width  <- mean(recalibrated_test$recal_upper -
                          recalibrated_test$recal_lower, na.rm = TRUE)
    recal_winkler <- mean(
      recal_width +
        (2 / alpha) * pmax(0, recalibrated_test$recal_lower - test_set$observed) +
        (2 / alpha) * pmax(0, test_set$observed - recalibrated_test$recal_upper),
      na.rm = TRUE
    )

    # ── Collect metrics ───────────────────────────────────────────────────
    metrics <- tibble(
      dataset     = bt$dataset,
      method      = c("Parametric", "Conformal", "Recalibrated"),
      coverage_95 = c(para_cov95, conf_cov95, recal_cov95),
      avg_width   = c(para_width, conf_width, recal_width),
      winkler     = c(para_winkler, conf_winkler, recal_winkler)
    )

    calibration_comparison[[nm]] <- list(
      dataset      = bt$dataset,
      metrics      = metrics,
      parametric   = parametric_test,
      conformal    = conformal_test,
      recalibrated = recalibrated_test,
      conf_radius  = conf_radius,
      recal_factor = recal_factor,
      n_cal        = n_cal,
      n_test       = n_origins - n_cal
    )

    cat(sprintf("      Parametric   95%% cov: %.1f%%, width: %.4f, Winkler: %.4f\n",
                para_cov95 * 100, para_width, para_winkler))
    cat(sprintf("      Conformal    95%% cov: %.1f%%, width: %.4f, Winkler: %.4f\n",
                conf_cov95 * 100, conf_width, conf_winkler))
    cat(sprintf("      Recalibrated 95%% cov: %.1f%%, width: %.4f, Winkler: %.4f\n",
                recal_cov95 * 100, recal_width, recal_winkler))
  },
  error = function(e) {
    cat(sprintf("    WARNING: Comparison failed for %s: %s\n", nm, e$message))
  })
}

###############################################################################
# Section D: Theoretical variance decomposition
###############################################################################

cat("  [D] Variance decomposition...\n")

# For each country, estimate V_captured / V_total
# V_captured: predicted variance from model's prediction interval width
# V_total: observed MSE from backtest residuals
# If ratio < 1: confirms underdispersion (model underestimates true variance)

variance_decomp <- list()

for (nm in names(benchmark$backtest_results)) {
  bt <- benchmark$backtest_results[[nm]]
  bt_tibble <- bt$result

  # Skip engines without prediction intervals
  has_intervals <- all(c("lower", "upper") %in% names(bt_tibble)) &&
    any(!is.na(bt_tibble$lower)) && any(!is.na(bt_tibble$upper))

  if (!has_intervals) {
    cat(sprintf("    %s/%s: no prediction intervals; skipping\n",
                bt$dataset, bt$engine))
    next
  }

  tryCatch({
    # Observed MSE per horizon
    mse_by_horizon <- bt_tibble |>
      filter(!is.na(lower), !is.na(upper), !is.na(predicted), !is.na(observed)) |>
      group_by(horizon) |>
      summarise(
        observed_mse  = mean((predicted - observed)^2, na.rm = TRUE),
        # Predicted variance: implied from 95% CI width
        # Under normality, 95% CI width = 2 * 1.96 * sigma
        predicted_var = mean(((upper - lower) / (2 * 1.96))^2, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        variance_ratio = predicted_var / observed_mse,
        dataset        = bt$dataset,
        engine         = bt$engine,
        # ratio < 1 indicates underdispersion
        underdispersed = variance_ratio < 1
      )

    variance_decomp[[nm]] <- mse_by_horizon

    cat(sprintf("    %s/%s: variance ratio = %.2f (h=7d), %.2f (h=28d) %s\n",
                bt$dataset, bt$engine,
                mse_by_horizon$variance_ratio[1],
                tail(mse_by_horizon$variance_ratio, 1),
                ifelse(all(mse_by_horizon$underdispersed), "[UNDERDISPERSED]", "")))
  },
  error = function(e) {
    cat(sprintf("    WARNING: Variance decomposition failed for %s: %s\n",
                nm, e$message))
  })
}

# Missing variance sources (documented):
# 1. Temporal correlation in δ(t): MLR assumes constant growth advantage δ,
#    but in practice δ varies over time (e.g., due to changing behavior, policy).
#    This time-varying δ introduces variance not captured by Fisher information.
# 2. Sampling correlation: genomic surveillance samples are not independent;
#    geographic and temporal clustering inflates variance.
# 3. Model misspecification: multinomial logistic regression assumes linearity
#    in log-odds, but true dynamics may be nonlinear (e.g., frequency-dependent
#    selection, heterogeneous contact patterns).

variance_decomp_df <- bind_rows(variance_decomp)

###############################################################################
# Save all results
###############################################################################

calibration_results <- list(
  pit_results            = pit_results,
  reliability            = reliability_df,
  calibration_comparison = calibration_comparison,
  variance_decomposition = variance_decomp_df
)

saveRDS(calibration_results, "analysis/results/calibration_comparison.rds")
cat("  Saved calibration_comparison.rds\n")

cat("[03_calibration] Complete.\n")
