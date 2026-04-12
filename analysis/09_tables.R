###############################################################################
# 09_tables.R — Nature Methods tables (LaTeX booktabs)
# lineagefreq validation analysis
###############################################################################

cat("[09_tables] Starting...\n")
source("analysis/00_setup.R")

benchmark   <- readRDS("analysis/results/benchmark_multicountry.rds")
calibration <- readRDS("analysis/results/calibration_comparison.rds")
decision    <- readRDS("analysis/results/decision_impact.rds")
sample_size <- tryCatch(readRDS("analysis/results/sample_size_analysis.rds"), error = function(e) NULL)

###############################################################################
# Table 1: Multi-country benchmark
###############################################################################

cat("  [Table 1] Multi-country benchmark...\n")

table1 <- benchmark$metrics |>
  filter(!is.na(mae)) |>
  select(dataset, engine, horizon, mae, median_ae, crps) |>
  mutate(
    mae       = sprintf("%.1f%%", mae * 100),
    median_ae = sprintf("%.1f%%", median_ae * 100),
    crps      = ifelse(is.na(crps), "—", sprintf("%.4f", crps))
  ) |>
  pivot_wider(
    names_from = horizon,
    values_from = c(mae, median_ae, crps),
    names_glue = "{.value}_{horizon}d"
  ) |>
  select(Dataset = dataset, Engine = engine,
         `MAE 7d` = mae_7d, `MAE 14d` = mae_14d, `MAE 28d` = mae_28d,
         `CRPS 7d` = crps_7d)

# Add Bedford Lab reference row
bedford_row <- tibble(
  Dataset    = "Bedford Lab*",
  Engine     = "mlr",
  `MAE 7d`   = "3.0%",
  `MAE 14d`  = "6.0%",
  `MAE 28d`  = "9.0%",
  `CRPS 7d`  = "—"
)

table1 <- bind_rows(table1, bedford_row)

save_table(table1, "table1",
           caption = "Multi-country forecasting benchmark. *Bedford Lab values from Abousamra et al.\\ 2024, PLOS Computational Biology.")

###############################################################################
# Table 2: Calibration comparison
###############################################################################

cat("  [Table 2] Calibration comparison...\n")

cal_metrics <- map_dfr(calibration$calibration_comparison, function(cc) {
  cc$metrics
}) |>
  filter(!is.na(coverage_95))

# Add KS D from PIT results
ks_data <- map_dfr(calibration$pit_results, function(pr) {
  tibble(dataset = pr$dataset, ks_D = pr$ks_D)
}) |>
  distinct(dataset, .keep_all = TRUE)

table2 <- cal_metrics |>
  left_join(ks_data, by = "dataset") |>
  mutate(
    coverage_95 = sprintf("%.1f%%", coverage_95 * 100),
    avg_width   = sprintf("%.4f", avg_width),
    winkler     = sprintf("%.4f", winkler),
    ks_D        = sprintf("%.3f", ks_D)
  ) |>
  select(Dataset = dataset, Method = method,
         `Coverage (95%)` = coverage_95,
         `Avg Width` = avg_width,
         `Winkler Score` = winkler,
         `KS D` = ks_D)

save_table(table2, "table2",
           caption = "Calibration comparison across prediction interval methods.")

###############################################################################
# Table 3: Decision impact
###############################################################################

cat("  [Table 3] Decision impact...\n")

if (nrow(decision) > 0) {
  table3 <- decision |>
    mutate(
      trigger_date        = format(trigger_date, "%Y-%m-%d"),
      trigger_date        = ifelse(is.na(trigger_date), "Not triggered", trigger_date),
      error_days          = ifelse(is.na(error_days), "—",
                                   sprintf("%+d", as.integer(error_days))),
      false_trigger_rate  = ifelse(is.na(false_trigger_rate), "—",
                                   sprintf("%.1f%%", false_trigger_rate * 100))
    ) |>
    select(Dataset = dataset, Method = method,
           `Trigger Date` = trigger_date,
           `Days vs Optimal` = error_days,
           `False Trigger Rate` = false_trigger_rate)

  save_table(table3, "table3",
             caption = "Decision impact: vaccine update trigger timing by method.")
} else {
  cat("    WARNING: No decision impact data; skipping Table 3\n")
}

###############################################################################
# Table 4: Sample size requirements
###############################################################################

cat("  [Table 4] Sample size requirements...\n")

if (!is.null(sample_size)) {
  table4 <- sample_size$summary |>
    mutate(
      mae_str     = sprintf("%.1f%% (±%.1f%%)", mae_mean * 100, mae_sd * 100),
      cov95_str   = sprintf("%.1f%% (±%.1f%%)", cov95_mean * 100, cov95_sd * 100),
      meets_mae   = ifelse(meets_mae5 > 0.5, "Yes", "No"),
      meets_cov   = ifelse(meets_cov90 > 0.5, "Yes", "No")
    ) |>
    select(
      `Seq/period`    = sample_size,
      `MAE`           = mae_str,
      `Coverage (95%)` = cov95_str,
      `MAE < 5%?`     = meets_mae,
      `Coverage > 90%?` = meets_cov
    )

  save_table(table4, "table4",
             caption = "Sample size requirements for forecast accuracy and calibration. 50 bootstrap replicates per sample size level.")
} else {
  cat("    WARNING: No sample size data; skipping Table 4\n")
}

cat("[09_tables] Complete.\n")
