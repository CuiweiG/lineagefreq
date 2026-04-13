###############################################################################
# Engine Calibration Comparison — Denmark BA.2
# Gap 2: Head-to-head calibration of all 5 engines
#
# Run from package root:
#   Rscript inst/extdata/engine_comparison/run_comparison.R
#
# Expected runtime: 2-6 hours for GARW (96-core EPYC)
# Intermediate results saved after each engine completes
###############################################################################

cat("═══════════════════════════════════════════════════════════\n")
cat("  Engine Calibration Comparison — Denmark BA.2\n")
cat(sprintf("  Started: %s\n", Sys.time()))
cat("═══════════════════════════════════════════════════════════\n\n")

library(lineagefreq)
library(dplyr)
library(tibble)

out_dir <- "inst/extdata/engine_comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load and prepare Denmark BA.2 data ────────────────────────────────────

ecdc <- readRDS("analysis/results/ecdc_prepared.rds")
dk_raw <- ecdc[["Denmark_BA2"]]

if (is.null(dk_raw)) stop("Denmark_BA2 not found in ecdc_prepared.rds")

# Collapse rare lineages: keep those with >5% peak frequency
dk <- collapse_lineages(dk_raw, min_freq = 0.05)

cat(sprintf("Denmark BA.2: %d dates, %d lineages after collapsing\n",
            length(unique(dk$.date)),
            length(attr(dk, "lineages"))))
cat("Lineages:", paste(attr(dk, "lineages"), collapse = ", "), "\n\n")

saveRDS(dk, file.path(out_dir, "denmark_ba2_collapsed.rds"))

# ── Backtest parameters ──────────────────────────────────────────────────

horizons  <- c(7L, 14L, 21L, 28L)
min_train <- 42L  # 6 weeks

# ── Engine configurations ────────────────────────────────────────────────

engines_config <- list(
  mlr = list(
    engine = "mlr",
    args = list()
  ),
  hier_mlr = list(
    engine = "hier_mlr",
    args = list()
  ),
  piantham = list(
    engine = "piantham",
    args = list(generation_time = 5.5)
  ),
  fga = list(
    engine = "fga",
    args = list(chains = 4, iter_warmup = 1000, iter_sampling = 1000)
  ),
  garw = list(
    engine = "garw",
    args = list(chains = 4, iter_warmup = 1000, iter_sampling = 1000)
  )
)

# ── Run backtests ────────────────────────────────────────────────────────

all_backtests <- list()
all_timings   <- list()

for (eng_name in names(engines_config)) {
  result_file <- file.path(out_dir, paste0("bt_", eng_name, ".rds"))

  # Skip if already computed
  if (file.exists(result_file)) {
    cat(sprintf("[%s] Loading cached result from %s\n", eng_name, result_file))
    all_backtests[[eng_name]] <- readRDS(result_file)
    all_timings[[eng_name]] <- list(engine = eng_name, elapsed = NA,
                                    status = "cached")
    next
  }

  cat(sprintf("[%s] Starting backtest at %s...\n", eng_name, Sys.time()))

  t0 <- proc.time()

  bt <- tryCatch({
    args <- c(
      list(data = dk, engines = eng_name,
           horizons = horizons, min_train = min_train),
      engines_config[[eng_name]]$args
    )
    do.call(backtest, args)
  },
  error = function(e) {
    cat(sprintf("[%s] FAILED: %s\n", eng_name, e$message))
    NULL
  })

  elapsed <- (proc.time() - t0)["elapsed"]

  if (!is.null(bt)) {
    all_backtests[[eng_name]] <- bt
    saveRDS(bt, result_file)
    cat(sprintf("[%s] Complete: %.1f sec, %d forecast-obs pairs\n",
                eng_name, elapsed, nrow(bt)))
  }

  all_timings[[eng_name]] <- list(
    engine  = eng_name,
    elapsed = elapsed,
    status  = if (!is.null(bt)) "success" else "failed",
    n_rows  = if (!is.null(bt)) nrow(bt) else 0
  )

  cat("\n")
}

# ── Compute calibration metrics for each engine ──────────────────────────

cat("\n═══════════════════════════════════════════════════════════\n")
cat("  Computing calibration metrics\n")
cat("═══════════════════════════════════════════════════════════\n\n")

cal_results <- list()
score_results <- list()

for (eng_name in names(all_backtests)) {
  bt <- all_backtests[[eng_name]]
  if (is.null(bt)) next

  # PIT / KS / reliability
  cal <- tryCatch(calibrate(bt), error = function(e) {
    cat(sprintf("[%s] calibrate() failed: %s\n", eng_name, e$message))
    NULL
  })

  if (!is.null(cal)) {
    cal_results[[eng_name]] <- list(
      engine     = eng_name,
      pit_values = cal$pit_values,
      ks_D       = cal$ks_test$statistic,
      ks_p       = cal$ks_test$p_value,
      n          = cal$n,
      reliability = cal$reliability
    )

    cat(sprintf("[%s] KS D = %.3f, p = %.2g, n = %d\n",
                eng_name, cal$ks_test$statistic, cal$ks_test$p_value, cal$n))
  }

  # Proper scoring rules
  scores <- tryCatch(score_forecasts(bt), error = function(e) {
    cat(sprintf("[%s] score_forecasts() failed: %s\n", eng_name, e$message))
    NULL
  })

  if (!is.null(scores)) {
    score_results[[eng_name]] <- scores
  }

  # Coverage at multiple levels
  bt_valid <- bt[!is.na(bt$lower) & !is.na(bt$upper), ]
  if (nrow(bt_valid) > 0) {
    sigma <- pmax((bt_valid$upper - bt_valid$lower) / (2 * qnorm(0.975)),
                   1e-10)
    z <- (bt_valid$observed - bt_valid$predicted) / sigma

    for (level in c(0.50, 0.80, 0.90, 0.95)) {
      z_crit <- qnorm(1 - (1 - level) / 2)
      cov <- mean(abs(z) <= z_crit, na.rm = TRUE)
      cat(sprintf("[%s] Coverage at %d%%: %.1f%%\n",
                  eng_name, level * 100, cov * 100))
    }
  }

  # Joint calibration (energy score)
  jcal <- tryCatch(calibrate_joint(bt), error = function(e) NULL)
  if (!is.null(jcal)) {
    cat(sprintf("[%s] Mean energy score: %.4f\n",
                eng_name, jcal$mean_energy_score))
    cal_results[[eng_name]]$energy_score <- jcal$mean_energy_score
  }

  cat("\n")
}

# ── Save consolidated results ────────────────────────────────────────────

cat("═══════════════════════════════════════════════════════════\n")
cat("  Saving consolidated results\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Summary table
summary_rows <- list()
for (eng_name in names(cal_results)) {
  cr <- cal_results[[eng_name]]
  sc <- score_results[[eng_name]]
  tm <- all_timings[[eng_name]]

  # Extract scores at 14d horizon
  mae_14 <- NA_real_
  crps_14 <- NA_real_
  coverage_14 <- NA_real_
  wis_14 <- NA_real_

  if (!is.null(sc)) {
    mae_row <- sc |> filter(horizon == 14, metric == "mae")
    if (nrow(mae_row) > 0) mae_14 <- mae_row$value[1]

    crps_row <- sc |> filter(horizon == 14, metric == "crps")
    if (nrow(crps_row) > 0) crps_14 <- crps_row$value[1]

    cov_row <- sc |> filter(horizon == 14, metric == "coverage")
    if (nrow(cov_row) > 0) coverage_14 <- cov_row$value[1]

    wis_row <- sc |> filter(horizon == 14, metric == "wis")
    if (nrow(wis_row) > 0) wis_14 <- wis_row$value[1]
  }

  summary_rows[[eng_name]] <- tibble(
    engine       = eng_name,
    ks_D         = cr$ks_D,
    ks_p         = cr$ks_p,
    n_pit        = cr$n,
    mae_14d      = mae_14,
    crps_14d     = crps_14,
    coverage_14d = coverage_14,
    wis_14d      = wis_14,
    energy_score = cr$energy_score %||% NA_real_,
    elapsed_s    = tm$elapsed %||% NA_real_
  )
}

summary_df <- bind_rows(summary_rows)

cat("Summary:\n")
print(summary_df)

saveRDS(list(
  backtests    = all_backtests,
  calibration  = cal_results,
  scores       = score_results,
  summary      = summary_df,
  timings      = all_timings,
  data_info    = list(
    dataset    = "Denmark_BA2",
    n_dates    = length(unique(dk$.date)),
    lineages   = attr(dk, "lineages"),
    date_range = range(dk$.date)
  )
), file.path(out_dir, "engine_comparison_results.rds"))

cat(sprintf("\nResults saved to %s/engine_comparison_results.rds\n", out_dir))
cat(sprintf("Completed: %s\n", Sys.time()))
