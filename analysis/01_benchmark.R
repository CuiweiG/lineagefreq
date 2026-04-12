#!/usr/bin/env Rscript
# ============================================================
# 01_benchmark.R — Forecast accuracy benchmark
# ============================================================
# Replicates and extends the evaluation protocol from
# Abousamra, Figgins & Bedford (2024, PLOS Comp Bio).
#
# Rolling-origin backtesting on CDC BA.2 and JN.1 datasets
# with all available engines.
#
# HONEST DISCLOSURE: piantham and hier_mlr are wrappers around
# mlr and do NOT require Stan. fga and garw require CmdStan
# and will be skipped if unavailable. The script tests each
# engine before including it.
# ============================================================

cat("[01_benchmark] Starting...\n")

if (!exists("theme_nature")) source("analysis/00_setup.R")

# ---- Load data ----
data(cdc_ba2_transition)
data(cdc_sarscov2_jn1)

x_ba2 <- lfq_data(cdc_ba2_transition, lineage = lineage,
                   date = date, count = count)
x_jn1 <- lfq_data(cdc_sarscov2_jn1, lineage = lineage,
                   date = date, count = count)

# ---- Detect available engines ----
# Test each engine on a small simulated dataset to see which work.
cat("[01_benchmark] Testing available engines...\n")
test_sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.9),
  n_timepoints = 10, total_per_tp = 200, seed = 999)

candidate_engines <- c("mlr", "hier_mlr", "piantham", "fga", "garw")
available_engines <- character(0)

for (eng in candidate_engines) {
  ok <- tryCatch({
    args <- list(data = test_sim, engine = eng)
    if (eng == "piantham") args$generation_time <- 5
    do.call(fit_model, args)
    TRUE
  }, error = function(e) FALSE)

  status <- if (ok) "available" else "UNAVAILABLE"
  cat(sprintf("  %-10s: %s\n", eng, status))
  if (ok) available_engines <- c(available_engines, eng)
}

cat(sprintf("[01_benchmark] Using engines: %s\n",
            paste(available_engines, collapse = ", ")))

if (length(available_engines) == 0) {
  stop("No engines available. Cannot proceed.")
}

# Note which engines are missing and why
missing_engines <- setdiff(candidate_engines, available_engines)
if (length(missing_engines) > 0) {
  cat(sprintf("[01_benchmark] LIMITATION: engines %s unavailable ",
              paste(missing_engines, collapse = ", ")))
  cat("(likely requires CmdStan installation).\n")
  cat("  Results below benchmark ONLY the available engines.\n")
  cat("  To include Bayesian engines, install CmdStan:\n")
  cat("  cmdstanr::install_cmdstan()\n\n")
}

# ---- Run backtests ----
horizons <- c(7L, 14L, 21L, 28L)

cat("[01_benchmark] Running backtest on BA.2 dataset...\n")
bt_ba2 <- backtest(x_ba2, engines = available_engines,
                   horizons = horizons, min_train = 42,
                   generation_time = 3.2)  # for piantham if available
cat(sprintf("  BA.2 backtest: %d prediction rows\n", nrow(bt_ba2)))

cat("[01_benchmark] Running backtest on JN.1 dataset...\n")
bt_jn1 <- backtest(x_jn1, engines = available_engines,
                   horizons = horizons, min_train = 42,
                   generation_time = 5.0)
cat(sprintf("  JN.1 backtest: %d prediction rows\n", nrow(bt_jn1)))

# ---- Score with all metrics ----
all_metrics <- c("mae", "rmse", "coverage", "wis",
                 "crps", "log_score", "dss", "calibration")

cat("[01_benchmark] Computing forecast scores...\n")
sc_ba2 <- score_forecasts(bt_ba2, metrics = all_metrics)
sc_jn1 <- score_forecasts(bt_jn1, metrics = all_metrics)

# ---- Compute detailed per-origin statistics ----
# We need per-origin MAE for bootstrap confidence intervals.
compute_detail <- function(bt, dataset_name) {
  bt_v <- bt[!is.na(bt$observed) & !is.na(bt$predicted), ]
  rows <- list()
  for (eng in unique(bt_v$engine)) {
    for (h in unique(bt_v$horizon)) {
      sub <- bt_v[bt_v$engine == eng & bt_v$horizon == h, ]
      if (nrow(sub) == 0L) next

      ae <- abs(sub$predicted - sub$observed)
      cov95 <- mean(sub$observed >= sub$lower &
                      sub$observed <= sub$upper, na.rm = TRUE)

      # Per-origin MAE for CI computation
      origin_mae <- sub |>
        group_by(origin_date) |>
        summarise(mae = mean(abs(predicted - observed)),
                  .groups = "drop")

      rows <- c(rows, list(tibble(
        dataset     = dataset_name,
        engine      = eng,
        horizon     = h,
        median_ae   = median(ae),
        mean_ae     = mean(ae),
        mae_se      = sd(origin_mae$mae) / sqrt(nrow(origin_mae)),
        mae_lower   = mean(ae) - 1.96 * sd(origin_mae$mae) / sqrt(nrow(origin_mae)),
        mae_upper   = mean(ae) + 1.96 * sd(origin_mae$mae) / sqrt(nrow(origin_mae)),
        coverage_95 = cov95,
        n_origins   = length(unique(sub$origin_date)),
        n_pairs     = nrow(sub)
      )))
    }
  }
  bind_rows(rows)
}

detail_ba2 <- compute_detail(bt_ba2, "BA.2")
detail_jn1 <- compute_detail(bt_jn1, "JN.1")
detail_all <- bind_rows(detail_ba2, detail_jn1)

# ---- Print summary ----
cat("\n[01_benchmark] Results summary:\n")
cat("  BA.2 transition:\n")
print(detail_ba2 |> select(engine, horizon, median_ae, mean_ae, coverage_95),
      n = 50)
cat("\n  JN.1 emergence:\n")
print(detail_jn1 |> select(engine, horizon, median_ae, mean_ae, coverage_95),
      n = 50)

# ---- Save results ----
saveRDS(list(
  bt_ba2 = bt_ba2, bt_jn1 = bt_jn1,
  sc_ba2 = sc_ba2, sc_jn1 = sc_jn1,
  detail = detail_all,
  available_engines = available_engines,
  missing_engines = missing_engines
), "analysis/results/benchmark_ba2.rds")

# Save JN.1 separately for downstream scripts that need only JN.1
saveRDS(list(bt = bt_jn1, scores = sc_jn1),
        "analysis/results/benchmark_jn1.rds")

cat("[01_benchmark] Complete.\n\n")
