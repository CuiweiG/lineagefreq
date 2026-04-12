#!/usr/bin/env Rscript
# ============================================================
# 02_calibration.R — Calibration diagnostics
# ============================================================
# This is the KEY section of the validation. It demonstrates that
# standard MLR prediction intervals are systematically too narrow
# — a finding invisible to MAE-based evaluation.
#
# We compute PIT histograms, reliability diagrams, and conformal
# prediction intervals to show the problem and its solution.
# ============================================================

cat("[02_calibration] Starting...\n")

if (!exists("theme_nature")) source("analysis/00_setup.R")

# ---- Load benchmark results ----
bench <- readRDS("analysis/results/benchmark_ba2.rds")
bt_ba2 <- bench$bt_ba2
bt_jn1 <- bench$bt_jn1

# ---- Calibration diagnostics ----
cat("[02_calibration] Computing PIT values for BA.2...\n")
cal_ba2 <- calibrate(bt_ba2)

cat("[02_calibration] Computing PIT values for JN.1...\n")
cal_jn1 <- calibrate(bt_jn1)

cat(sprintf("  BA.2: KS D = %.4f, p = %s, n = %d\n",
  cal_ba2$ks_test$statistic,
  format.pval(cal_ba2$ks_test$p_value, digits = 3),
  cal_ba2$n))
cat(sprintf("  JN.1: KS D = %.4f, p = %s, n = %d\n",
  cal_jn1$ks_test$statistic,
  format.pval(cal_jn1$ks_test$p_value, digits = 3),
  cal_jn1$n))

# ---- Mean calibration error ----
mce_ba2 <- mean(abs(cal_ba2$reliability$observed -
                     cal_ba2$reliability$nominal))
mce_jn1 <- mean(abs(cal_jn1$reliability$observed -
                     cal_jn1$reliability$nominal))
cat(sprintf("  BA.2 mean calibration error: %.4f\n", mce_ba2))
cat(sprintf("  JN.1 mean calibration error: %.4f\n", mce_jn1))

# ---- Conformal prediction comparison ----
cat("[02_calibration] Running conformal prediction on JN.1...\n")
data(cdc_sarscov2_jn1)
x_jn1 <- lfq_data(cdc_sarscov2_jn1, lineage = lineage,
                   date = date, count = count)
fit_jn1 <- fit_model(x_jn1, engine = "mlr")

fc_parametric <- forecast(fit_jn1, horizon = 28, ci_level = 0.95)
fc_conformal  <- conformal_forecast(fit_jn1, x_jn1, horizon = 28,
                                     ci_level = 0.95, seed = 42)

# Compare interval widths
param_fc <- fc_parametric[fc_parametric$.type == "forecast", ]
conf_fc  <- fc_conformal[fc_conformal$.type == "forecast", ]
cat(sprintf("  Parametric interval width: %.4f\n",
            mean(param_fc$.upper - param_fc$.lower)))
cat(sprintf("  Conformal interval width:  %.4f\n",
            mean(conf_fc$.upper - conf_fc$.lower)))

# Also run conformal on BA.2
data(cdc_ba2_transition)
x_ba2 <- lfq_data(cdc_ba2_transition, lineage = lineage,
                   date = date, count = count)
fit_ba2 <- fit_model(x_ba2, engine = "mlr")
fc_conf_ba2 <- conformal_forecast(fit_ba2, x_ba2, horizon = 28,
                                   ci_level = 0.95, seed = 42)

# ---- Save results ----
saveRDS(list(
  cal_ba2 = cal_ba2,
  cal_jn1 = cal_jn1,
  mce_ba2 = mce_ba2,
  mce_jn1 = mce_jn1,
  fc_parametric = fc_parametric,
  fc_conformal  = fc_conformal,
  fc_conf_ba2   = fc_conf_ba2,
  fit_jn1 = fit_jn1,
  fit_ba2 = fit_ba2
), "analysis/results/calibration_results.rds")

cat("[02_calibration] Complete.\n\n")
