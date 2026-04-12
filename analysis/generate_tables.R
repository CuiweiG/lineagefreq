#!/usr/bin/env Rscript
# ============================================================
# Generate LaTeX tables from validation results
# ============================================================
#
# Reads analysis/results/*.rds and produces .tex files in
# analysis/tables/. Requires kableExtra.
#
# Usage:
#   source("analysis/generate_tables.R")
# ============================================================

library(kableExtra)
library(dplyr)
library(tibble)

cat("Generating LaTeX tables...\n")

# ---- Table 1: Benchmark comparison ----
bench <- readRDS("analysis/results/benchmark_raw.rds")
detail <- bench$detail

# Add Bedford Lab reference rows
bedford <- tibble(
  dataset = "BA.2 (Bedford Lab 2024)",
  engine = "mlr (reference)",
  horizon = c(7, 14, 21, 28),
  median_ae = c(0.006, NA, NA, NA),
  mean_ae = c(NA, 0.06, 0.08, 0.09),
  coverage_95 = rep(NA_real_, 4),
  n_origins = rep(NA_integer_, 4)
)

# Combine CRPS from scores
crps_ba2 <- bench$sc_ba2 |>
  filter(metric == "crps") |>
  select(engine, horizon, crps = value)

tab1_data <- detail |>
  filter(dataset == "BA.2") |>
  left_join(crps_ba2, by = c("engine", "horizon")) |>
  bind_rows(bedford) |>
  mutate(
    median_ae = sprintf("%.4f", median_ae),
    mean_ae = sprintf("%.4f", mean_ae),
    crps = ifelse(is.na(crps), "---", sprintf("%.4f", crps)),
    coverage_95 = ifelse(is.na(coverage_95), "---",
                         sprintf("%.1f%%", coverage_95 * 100)),
    n_origins = ifelse(is.na(n_origins), "---",
                       as.character(n_origins))
  ) |>
  select(Engine = engine, Horizon = horizon,
         `Median AE` = median_ae, `Mean AE` = mean_ae,
         CRPS = crps, `Coverage (95%)` = coverage_95,
         Origins = n_origins)

kbl(tab1_data, format = "latex", booktabs = TRUE,
    caption = "Forecast accuracy: lineagefreq vs Bedford Lab (2024)",
    label = "tab:benchmark") |>
  kable_styling(latex_options = c("hold_position")) |>
  save_kable("analysis/tables/table1_benchmark.tex")

cat("  Table 1 saved.\n")


# ---- Table 2: Calibration results ----
cal <- readRDS("analysis/results/calibration_raw.rds")

tab2_data <- tibble(
  Dataset = c("BA.2", "JN.1"),
  `KS statistic` = sprintf("%.4f",
    c(cal$cal_ba2$ks_test$statistic, cal$cal_jn1$ks_test$statistic)),
  `KS p-value` = format.pval(
    c(cal$cal_ba2$ks_test$p_value, cal$cal_jn1$ks_test$p_value),
    digits = 3),
  `Mean cal. error` = sprintf("%.4f",
    c(mean(abs(cal$cal_ba2$reliability$observed -
               cal$cal_ba2$reliability$nominal)),
      mean(abs(cal$cal_jn1$reliability$observed -
               cal$cal_jn1$reliability$nominal)))),
  `Coverage at 90%` = sprintf("%.1f%%",
    c(cal$cal_ba2$reliability$observed[
        cal$cal_ba2$reliability$nominal == 0.9] * 100,
      cal$cal_jn1$reliability$observed[
        cal$cal_jn1$reliability$nominal == 0.9] * 100)),
  N = c(cal$cal_ba2$n, cal$cal_jn1$n)
)

kbl(tab2_data, format = "latex", booktabs = TRUE,
    caption = "Calibration diagnostics: PIT uniformity test",
    label = "tab:calibration") |>
  kable_styling(latex_options = c("hold_position")) |>
  save_kable("analysis/tables/table2_calibration.tex")

cat("  Table 2 saved.\n")


# ---- Table 3: Surveillance optimization ----
surv <- readRDS("analysis/results/surveillance_raw.rds")

tab3_data <- tibble(
  Strategy = c("Thompson sampling", "UCB"),
  `Total sequences` = c(
    sum(surv$ad_thompson$allocations$n_allocated),
    sum(surv$ad_ucb$allocations$n_allocated)),
  `Mean uncertainty` = sprintf("%.6f", c(
    mean(surv$ad_thompson$summary$mean_uncertainty),
    mean(surv$ad_ucb$summary$mean_uncertainty))),
  Rounds = c(
    length(unique(surv$ad_thompson$allocations$round)),
    length(unique(surv$ad_ucb$allocations$round)))
)

kbl(tab3_data, format = "latex", booktabs = TRUE,
    caption = "Adaptive allocation strategy comparison",
    label = "tab:surveillance") |>
  kable_styling(latex_options = c("hold_position")) |>
  save_kable("analysis/tables/table3_surveillance.tex")

cat("  Table 3 saved.\n")
cat("All tables generated.\n")
