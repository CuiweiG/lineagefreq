###############################################################################
# 11_multitransition_validation.R
# Comprehensive multi-transition calibration assessment
# Validates lineagefreq + CAPS across 12 US transitions (2021-2026)
#
# Run from package root:
#   source("analysis/11_multitransition_validation.R")
#
# Requires: analysis/data/cdc_variant_proportions_full.csv (379 MB)
# Expected runtime: 15-30 minutes on 96-core EPYC
###############################################################################

cat("[11_multitransition] Starting...\n")
source("analysis/00_setup.R")
library(lineagefreq)
library(data.table)

set.seed(20260413)

# ─── Load CDC data ────────────────────────────────────────────────────────

cat("  Loading CDC variant proportions...\n")
cdc <- fread("analysis/data/cdc_variant_proportions_full.csv")
cat(sprintf("  %d rows, columns: %s\n", nrow(cdc),
            paste(head(names(cdc), 8), collapse = ", ")))

# Parse date — week_ending may be character or POSIXct
cdc[, date := as.Date(week_ending)]
cat(sprintf("  Date range: %s to %s\n", min(cdc$date), max(cdc$date)))

# Filter to USA national, weighted model, biweekly
cdc_usa <- cdc[usa_or_hhsregion == "USA"]
if ("modeltype" %in% names(cdc_usa)) {
  cdc_usa <- cdc_usa[modeltype == "weighted" | is.na(modeltype)]
}
cat(sprintf("  National US data: %d rows\n", nrow(cdc_usa)))

# ─── Define 12 major US transitions ──────────────────────────────────────

transitions <- list(
  list(name = "T01_Alpha_rise",
       start = "2021-01-01", end = "2021-06-01",
       desc = "B.1.1.7 (Alpha) replacing pre-Alpha lineages"),
  list(name = "T02_Alpha_to_Delta",
       start = "2021-04-01", end = "2021-09-01",
       desc = "B.1.617.2 (Delta) replacing Alpha"),
  list(name = "T03_Delta_to_Omicron",
       start = "2021-10-01", end = "2022-02-01",
       desc = "BA.1 (Omicron) replacing Delta"),
  list(name = "T04_BA1_to_BA2",
       start = "2021-12-01", end = "2022-06-30",
       desc = "BA.2 replacing BA.1 (previously validated)"),
  list(name = "T05_BA2_to_BA45",
       start = "2022-04-01", end = "2022-10-01",
       desc = "BA.4/5 replacing BA.2/BA.2.12.1"),
  list(name = "T06_BA5_to_BQ_XBB",
       start = "2022-09-01", end = "2023-04-01",
       desc = "BQ.1.1 and XBB.1.5 replacing BA.5"),
  list(name = "T07_XBB15_dominance",
       start = "2023-02-01", end = "2023-09-01",
       desc = "XBB.1.5 dominance with EG.5 emergence"),
  list(name = "T08_EG5_to_JN1",
       start = "2023-08-01", end = "2024-03-01",
       desc = "JN.1 replacing EG.5/HV.1 (previously validated)"),
  list(name = "T09_JN1_to_KP",
       start = "2024-02-01", end = "2024-09-01",
       desc = "KP.2/KP.3 replacing JN.1"),
  list(name = "T10_KP3_to_XEC",
       start = "2024-07-01", end = "2025-02-01",
       desc = "KP.3.1.1 and XEC replacing KP.3"),
  list(name = "T11_XEC_to_XFG",
       start = "2025-01-01", end = "2025-09-01",
       desc = "XFG replacing XEC/LP.8.1"),
  list(name = "T12_XFG_sublineages",
       start = "2025-08-01", end = "2026-04-01",
       desc = "XFG.1.1/NB.1.8.1 diversification (most recent)")
)

cat(sprintf("  Processing %d transitions...\n", length(transitions)))

# ─── Assumed total sequences per biweek ──────────────────────────────────
# CDC reports proportions, not counts. We reconstruct counts assuming
# ~8000 sequences per national biweekly period (conservative estimate).
ASSUMED_TOTAL <- 8000L

# ─── Process each transition ─────────────────────────────────────────────

transition_results <- list()

for (tr in transitions) {
  cat(sprintf("\n  -- %s: %s --\n", tr$name, tr$desc))

  # Filter time window
  tr_data <- cdc_usa[date >= as.Date(tr$start) & date <= as.Date(tr$end)]

  if (nrow(tr_data) == 0) {
    cat("    No data in this window, skipping\n")
    next
  }

  # Keep variants with >5% peak share
  variant_peaks <- tr_data[, .(peak = max(share, na.rm = TRUE)), by = variant]
  major_variants <- variant_peaks[peak >= 0.05]$variant
  tr_data <- tr_data[variant %in% major_variants]

  n_dates    <- length(unique(tr_data$date))
  n_lineages <- length(major_variants)

  cat(sprintf("    %d rows, %d variants, %d time points\n",
              nrow(tr_data), n_lineages, n_dates))

  if (n_dates < 6 || n_lineages < 2) {
    cat("    Insufficient data, skipping\n")
    next
  }

  # Reconstruct counts from proportions
  tr_counts <- tr_data[, .(
    lineage = variant,
    date    = date,
    count   = as.integer(round(share * ASSUMED_TOTAL))
  )]
  tr_counts <- tr_counts[count > 0]

  # Create lfq_data
  lfq <- tryCatch(
    lfq_data(tr_counts, lineage = lineage, date = date, count = count),
    error = function(e) {
      cat(sprintf("    lfq_data failed: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(lfq)) next

  # Fit MLR
  fit <- tryCatch(
    fit_model(lfq, engine = "mlr"),
    error = function(e) {
      cat(sprintf("    fit_model failed: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(fit)) next

  # Backtest — returns lfq_backtest tibble directly
  bt <- tryCatch(
    backtest(lfq, engines = "mlr",
             horizons = c(7L, 14L, 21L, 28L), min_train = 42L),
    error = function(e) {
      cat(sprintf("    backtest failed: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(bt)) next

  cat(sprintf("    Backtest: %d rows, %d origins\n",
              nrow(bt), attr(bt, "n_origins")))

  # Calibration
  cal <- tryCatch(calibrate(bt), error = function(e) NULL)

  ks_D <- NA_real_; ks_p <- NA_real_
  if (!is.null(cal)) {
    ks_D <- cal$ks_test$statistic
    ks_p <- cal$ks_test$p_value
  }

  # Coverage by horizon
  bt_valid <- bt[!is.na(bt$lower) & !is.na(bt$upper) & !is.na(bt$observed), ]
  coverage_by_h <- bt_valid |>
    dplyr::group_by(horizon) |>
    dplyr::summarise(
      coverage = mean(observed >= lower & observed <= upper, na.rm = TRUE),
      mae      = mean(abs(predicted - observed), na.rm = TRUE),
      n        = dplyr::n(),
      .groups  = "drop"
    )

  # Variance ratios by horizon
  vr_by_h <- bt_valid |>
    dplyr::group_by(horizon) |>
    dplyr::summarise(
      var_ratio = mean(((upper - lower) / (2 * 1.96))^2, na.rm = TRUE) /
                  mean((predicted - observed)^2, na.rm = TRUE),
      .groups = "drop"
    )

  # CAPS — pass fit + data, not calibration_data
  caps <- tryCatch(
    caps_forecast(fit, lfq,
                  horizons = c(7L, 14L, 21L, 28L),
                  method = "caps", alpha = 0.05),
    error = function(e) {
      cat(sprintf("    CAPS failed: %s\n", conditionMessage(e)))
      NULL
    }
  )

  # CAPS coverage via evaluate_caps
  caps_eval <- NULL
  if (!is.null(caps)) {
    caps_eval <- tryCatch(
      evaluate_caps(caps, bt),
      error = function(e) {
        cat(sprintf("    evaluate_caps failed: %s\n", conditionMessage(e)))
        NULL
      }
    )
  }

  avg_coverage <- mean(coverage_by_h$coverage, na.rm = TRUE)
  avg_mae      <- mean(coverage_by_h$mae, na.rm = TRUE)
  avg_vr       <- mean(vr_by_h$var_ratio, na.rm = TRUE)
  avg_caps_cov <- if (!is.null(caps_eval)) mean(caps_eval$caps_cov, na.rm = TRUE) else NA_real_

  cat(sprintf("    KS D=%.3f | Avg 95%% cov=%.1f%% | CAPS cov=%.1f%% | MAE=%.1f pp | R=%.2f\n",
              ks_D, avg_coverage * 100,
              if (is.na(avg_caps_cov)) NA else avg_caps_cov * 100,
              avg_mae * 100, avg_vr))

  transition_results[[tr$name]] <- list(
    name        = tr$name,
    description = tr$desc,
    period      = paste(tr$start, "to", tr$end),
    n_variants  = n_lineages,
    n_dates     = n_dates,
    n_origins   = attr(bt, "n_origins"),
    ks_D        = ks_D,
    ks_p        = ks_p,
    coverage_by_horizon = coverage_by_h,
    variance_ratios     = vr_by_h,
    caps_eval           = caps_eval,
    avg_coverage        = avg_coverage,
    avg_mae             = avg_mae,
    avg_var_ratio       = avg_vr,
    avg_caps_coverage   = avg_caps_cov,
    caps_R_hat          = if (!is.null(caps)) caps$R_hat else NULL,
    caps_psi_K          = if (!is.null(caps)) caps$psi_K else NULL
  )
}

# ─── Summary ─────────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 60), "\n")
cat("Multi-transition calibration summary (US national, 2021-2026)\n")
cat(strrep("=", 60), "\n\n")

summary_rows <- lapply(transition_results, function(r) {
  data.frame(
    transition = r$name,
    period     = r$period,
    variants   = r$n_variants,
    origins    = r$n_origins,
    KS_D       = round(r$ks_D, 3),
    para_cov   = sprintf("%.1f%%", r$avg_coverage * 100),
    caps_cov   = sprintf("%.1f%%", r$avg_caps_coverage * 100),
    MAE_pp     = sprintf("%.1f", r$avg_mae * 100),
    avg_R      = sprintf("%.2f", r$avg_var_ratio),
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
print(summary_df, row.names = FALSE)

# Universal miscalibration test
all_ks  <- sapply(transition_results, function(r) r$ks_D)
all_cov <- sapply(transition_results, function(r) r$avg_coverage)
all_caps <- sapply(transition_results, function(r) r$avg_caps_coverage)

cat("\n-- Universal miscalibration test --\n")
cat(sprintf("Transitions evaluated: %d\n", length(all_ks)))
cat(sprintf("All KS D > 0.1: %s\n",
            ifelse(all(all_ks > 0.1, na.rm = TRUE), "YES", "NO")))
cat(sprintf("All parametric coverage < 95%%: %s\n",
            ifelse(all(all_cov < 0.95, na.rm = TRUE), "YES", "NO")))
cat(sprintf("KS D range: %.3f to %.3f\n",
            min(all_ks, na.rm = TRUE), max(all_ks, na.rm = TRUE)))
cat(sprintf("Parametric coverage range: %.1f%% to %.1f%%\n",
            min(all_cov, na.rm = TRUE) * 100, max(all_cov, na.rm = TRUE) * 100))
cat(sprintf("CAPS coverage range: %.1f%% to %.1f%%\n",
            min(all_caps, na.rm = TRUE) * 100, max(all_caps, na.rm = TRUE) * 100))

saveRDS(transition_results, "analysis/results/multitransition_results.rds")
saveRDS(summary_df, "analysis/results/multitransition_summary.rds")

cat("\n[11_multitransition] Complete.\n")
