#!/usr/bin/env Rscript
# ============================================================
# 03_surveillance.R — Surveillance optimization analysis
# ============================================================
# Demonstrates adaptive allocation, EVOI, detection horizon,
# and sequential alerting. Documents SPRT limitations on
# biweekly data honestly.
# ============================================================

cat("[03_surveillance] Starting...\n")

if (!exists("theme_nature")) source("analysis/00_setup.R")

# ---- EVOI analysis ----
cat("[03_surveillance] Computing Expected Value of Information...\n")
data(cdc_ba2_transition)
x_ba2 <- lfq_data(cdc_ba2_transition, lineage = lineage,
                   date = date, count = count)
fit_ba2 <- fit_model(x_ba2, engine = "mlr")

ev <- surveillance_value(fit_ba2, n_current = 500,
                         n_additional = seq(10L, 500L, by = 10L))
cat(sprintf("  Current variance: %.6f\n", ev$current_uncertainty))

# ---- Adaptive allocation ----
cat("[03_surveillance] Running adaptive allocation (Thompson)...\n")
ad_thompson <- adaptive_design(x_ba2, capacity = 200,
                                n_rounds = 10,
                                strategy = "thompson", seed = 42)

cat("[03_surveillance] Running adaptive allocation (UCB)...\n")
ad_ucb <- adaptive_design(x_ba2, capacity = 200,
                           n_rounds = 10,
                           strategy = "ucb", seed = 42)

# ---- Detection horizon ----
cat("[03_surveillance] Computing detection horizons...\n")
dh_fast <- detection_horizon(initial_prev = 0.001, growth_rate = 1.5,
                              n_per_period = 500, n_periods = 20)
dh_slow <- detection_horizon(initial_prev = 0.001, growth_rate = 1.1,
                              n_per_period = 500, n_periods = 30)
dh_low  <- detection_horizon(initial_prev = 0.001, growth_rate = 1.3,
                              n_per_period = 50, n_periods = 30)

cat(sprintf("  Fast growth (1.5x/wk, 500 seq): %s weeks to detect\n",
  ifelse(is.na(attr(dh_fast, "weeks_to_detection")), "not reached",
         as.character(attr(dh_fast, "weeks_to_detection")))))
cat(sprintf("  Slow growth (1.1x/wk, 500 seq): %s weeks to detect\n",
  ifelse(is.na(attr(dh_slow, "weeks_to_detection")), "not reached",
         as.character(attr(dh_slow, "weeks_to_detection")))))
cat(sprintf("  Low resource (1.3x/wk, 50 seq): %s weeks to detect\n",
  ifelse(is.na(attr(dh_low, "weeks_to_detection")), "not reached",
         as.character(attr(dh_low, "weeks_to_detection")))))

# ---- Alert threshold on JN.1 ----
cat("[03_surveillance] Running SPRT on JN.1 data...\n")
data(cdc_sarscov2_jn1)
x_jn1 <- lfq_data(cdc_sarscov2_jn1, lineage = lineage,
                   date = date, count = count)

alerts_sprt  <- alert_threshold(x_jn1, method = "sprt",
                                 delta_1 = 0.02, alpha = 0.05)
alerts_cusum <- alert_threshold(x_jn1, method = "cusum",
                                 threshold = 3)

jn1_alert <- alerts_sprt[alerts_sprt$lineage == "JN.1", ]
cat(sprintf("  JN.1 SPRT alert: %s (stat = %.4f)\n",
  jn1_alert$alert, jn1_alert$statistic))

# When did JN.1 cross 5%?
jn1_freq <- cdc_sarscov2_jn1[cdc_sarscov2_jn1$lineage == "JN.1", ]
cross_5pct <- min(jn1_freq$date[jn1_freq$proportion >= 0.05])
cat(sprintf("  JN.1 crossed 5%%: %s\n", cross_5pct))

# HONEST DOCUMENTATION: SPRT may not trigger on biweekly data.
# The test accumulates evidence one observation at a time, and
# with only 2-3 observations between JN.1's first detection and
# its 5% crossing, there is insufficient evidence to cross the
# SPRT boundary. This is a real limitation of sequential testing
# on low-frequency surveillance data.
if (!jn1_alert$alert) {
  cat("  NOTE: SPRT did not trigger on biweekly data.\n")
  cat("  This is expected: sequential tests need many observations.\n")
  cat("  Recommendation: use summarize_emerging() for biweekly data,\n")
  cat("  reserve alert_threshold() for weekly or more frequent data.\n")
}

# ---- Save results ----
saveRDS(list(
  ev = ev,
  ad_thompson = ad_thompson,
  ad_ucb = ad_ucb,
  dh_fast = dh_fast,
  dh_slow = dh_slow,
  dh_low = dh_low,
  alerts_sprt = alerts_sprt,
  alerts_cusum = alerts_cusum,
  cross_5pct = cross_5pct,
  sprt_triggered = jn1_alert$alert
), "analysis/results/surveillance_results.rds")

cat("[03_surveillance] Complete.\n\n")
