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

  tryCatch({
    cat(sprintf("    Calibrating %s/%s...\n", bt$dataset, bt$engine))

    cal <- calibrate(bt$result)

    # Extract PIT values
    pit_values <- cal$pit_values

    # KS test for uniformity
    ks_result <- ks.test(pit_values, "punif")

    pit_results[[nm]] <- list(
      dataset    = bt$dataset,
      engine     = bt$engine,
      pit_values = pit_values,
      ks_D       = ks_result$statistic,
      ks_p       = ks_result$p.value,
      n_pit      = length(pit_values),
      calibration_obj = cal
    )

    cat(sprintf("      KS D = %.3f, p = %.3g (n = %d)\n",
                ks_result$statistic, ks_result$p.value, length(pit_values)))
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

nominal_levels <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

reliability_data <- list()

for (nm in names(benchmark$backtest_results)) {
  bt <- benchmark$backtest_results[[nm]]

  tryCatch({
    forecasts <- bt$result$forecasts

    obs_coverage <- sapply(nominal_levels, function(level) {
      col_lo <- paste0("lower_", level * 100)
      col_hi <- paste0("upper_", level * 100)

      if (all(c(col_lo, col_hi) %in% names(forecasts))) {
        mean(forecasts$observed >= forecasts[[col_lo]] &
             forecasts$observed <= forecasts[[col_hi]], na.rm = TRUE)
      } else {
        NA_real_
      }
    })

    reliability_data[[nm]] <- tibble(
      dataset          = bt$dataset,
      engine           = bt$engine,
      nominal_coverage = nominal_levels,
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

all_datasets <- c(
  benchmark$backtest_results |>
    Filter(function(x) x$engine == "mlr", x = _)
)

calibration_comparison <- list()

for (nm in names(all_datasets)) {
  bt <- all_datasets[[nm]]

  cat(sprintf("    Comparing methods for %s...\n", bt$dataset))

  tryCatch({
    data_obj <- bt$result

    # 1. Parametric forecast (default from backtest)
    parametric <- data_obj$forecasts |>
      mutate(method = "Parametric")

    # 2. Conformal forecast
    conformal_result <- tryCatch({
      conformal_forecast(data_obj)
    }, error = function(e) {
      cat(sprintf("      WARNING: conformal_forecast() failed: %s\n", e$message))
      NULL
    })

    # 3. Recalibrated forecast (isotonic regression)
    recalibrated_result <- tryCatch({
      recalibrate(data_obj)
    }, error = function(e) {
      cat(sprintf("      WARNING: recalibrate() failed: %s\n", e$message))
      NULL
    })

    # Compute metrics for each method
    compute_method_metrics <- function(fcast, method_name) {
      if (is.null(fcast)) return(NULL)

      fcast_df <- if (is.data.frame(fcast)) fcast else fcast$forecasts

      # 95% coverage
      cov95 <- tryCatch({
        mean(fcast_df$observed >= fcast_df$lower_95 &
             fcast_df$observed <= fcast_df$upper_95, na.rm = TRUE)
      }, error = function(e) NA_real_)

      # Average interval width at 95%
      avg_width <- tryCatch({
        mean(fcast_df$upper_95 - fcast_df$lower_95, na.rm = TRUE)
      }, error = function(e) NA_real_)

      # Winkler score (interval score) at 95%
      # Winkler (1972): S = width + (2/α)(lower - obs) if obs < lower
      #                         + (2/α)(obs - upper) if obs > upper
      alpha <- 0.05
      winkler <- tryCatch({
        width  <- fcast_df$upper_95 - fcast_df$lower_95
        below  <- pmax(0, fcast_df$lower_95 - fcast_df$observed)
        above  <- pmax(0, fcast_df$observed - fcast_df$upper_95)
        mean(width + (2 / alpha) * (below + above), na.rm = TRUE)
      }, error = function(e) NA_real_)

      tibble(
        dataset  = bt$dataset,
        method   = method_name,
        coverage_95 = cov95,
        avg_width   = avg_width,
        winkler     = winkler
      )
    }

    m_para   <- compute_method_metrics(parametric, "Parametric")
    m_conf   <- compute_method_metrics(conformal_result, "Conformal")
    m_recal  <- compute_method_metrics(recalibrated_result, "Recalibrated")

    calibration_comparison[[nm]] <- list(
      dataset      = bt$dataset,
      metrics      = bind_rows(m_para, m_conf, m_recal),
      parametric   = parametric,
      conformal    = conformal_result,
      recalibrated = recalibrated_result
    )

    cat(sprintf("      Parametric  95%% cov: %.1f%%, width: %.4f\n",
                m_para$coverage_95 * 100, m_para$avg_width))
    if (!is.null(m_conf))
      cat(sprintf("      Conformal   95%% cov: %.1f%%, width: %.4f\n",
                  m_conf$coverage_95 * 100, m_conf$avg_width))
    if (!is.null(m_recal))
      cat(sprintf("      Recalibrated 95%% cov: %.1f%%, width: %.4f\n",
                  m_recal$coverage_95 * 100, m_recal$avg_width))
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
# V_captured: predicted variance from model's Fisher information
# V_total: observed MSE from backtest residuals
# If ratio < 1: confirms underdispersion (model underestimates true variance)

variance_decomp <- list()

for (nm in names(benchmark$backtest_results)) {
  bt <- benchmark$backtest_results[[nm]]

  tryCatch({
    fcast <- bt$result$forecasts

    # Observed MSE per horizon
    mse_by_horizon <- fcast |>
      group_by(horizon) |>
      summarise(
        observed_mse   = mean((predicted - observed)^2, na.rm = TRUE),
        # Predicted variance: average of (upper_95 - lower_95)/(2*1.96))^2
        # Under normality, 95% CI width = 2 * 1.96 * sigma
        predicted_var  = mean(((upper_95 - lower_95) / (2 * 1.96))^2, na.rm = TRUE),
        variance_ratio = predicted_var / observed_mse,
        .groups = "drop"
      ) |>
      mutate(
        dataset = bt$dataset,
        engine  = bt$engine,
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
