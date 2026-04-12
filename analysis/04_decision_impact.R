###############################################################################
# 04_decision_impact.R — Decision impact analysis (NEW — Nature Methods)
# lineagefreq validation analysis
###############################################################################

cat("[04_decision_impact] Starting...\n")
source("analysis/00_setup.R")

set.seed(20240414)

benchmark  <- readRDS("analysis/results/benchmark_multicountry.rds")
calibration <- readRDS("analysis/results/calibration_comparison.rds")

###############################################################################
# Section A: Vaccine update trigger analysis
###############################################################################

cat("  [A] Vaccine update trigger analysis...\n")

# Scenario: A public health agency triggers a vaccine composition update review
# when a new variant's 95% LOWER confidence bound exceeds 30% frequency.
# This threshold is relevant because it signals the new variant is likely to
# become dominant, warranting urgent vaccine reformulation.
#
# We compare trigger timing across three interval estimation methods:
# - Parametric: standard MLR-based intervals (potentially too narrow)
# - Conformal: adaptive conformal intervals (wider, coverage-guaranteed)
# - Recalibrated: isotonic regression recalibrated intervals

TRIGGER_THRESHOLD <- 0.30  # 30% frequency

# Find BA.2 backtest results (US built-in data preferred)
ba2_keys <- names(benchmark$backtest_results) |>
  grep("BA2|ba2", x = _, value = TRUE)

if (length(ba2_keys) == 0) {
  cat("    WARNING: No BA.2 backtest results found; using first available\n")
  ba2_keys <- names(benchmark$backtest_results)[1]
}

trigger_results <- list()

for (key in ba2_keys) {
  bt <- benchmark$backtest_results[[key]]

  tryCatch({
    cat(sprintf("    Analyzing trigger timing for %s...\n", bt$dataset))

    fcast <- bt$result$forecasts

    # Identify BA.2 lineage column (may be named differently)
    ba2_lineage <- fcast |>
      filter(grepl("BA\\.2|BA2", lineage, ignore.case = TRUE)) |>
      pull(lineage) |>
      unique()

    if (length(ba2_lineage) == 0) {
      cat("      No BA.2 lineage found; using lineage with largest max frequency\n")
      ba2_lineage <- fcast |>
        group_by(lineage) |>
        summarise(max_obs = max(observed, na.rm = TRUE)) |>
        slice_max(max_obs, n = 1) |>
        pull(lineage)
    }

    ba2_fcast <- fcast |>
      filter(lineage == ba2_lineage[1])

    # Actual crossing date: first date when BA.2 observed > 30%
    actual_crossing <- ba2_fcast |>
      filter(observed > TRIGGER_THRESHOLD) |>
      slice_min(origin_date) |>
      pull(origin_date) |>
      min(na.rm = TRUE)

    # Parametric trigger: first origin where lower_95 > threshold
    parametric_trigger <- ba2_fcast |>
      filter(lower_95 > TRIGGER_THRESHOLD) |>
      slice_min(origin_date) |>
      pull(origin_date) |>
      min(na.rm = TRUE)

    # Conformal trigger
    conformal_cal <- calibration$calibration_comparison[[key]]
    conformal_trigger <- NA
    if (!is.null(conformal_cal$conformal)) {
      conf_fcast <- if (is.data.frame(conformal_cal$conformal)) {
        conformal_cal$conformal
      } else {
        conformal_cal$conformal$forecasts
      }

      if (!is.null(conf_fcast)) {
        conf_ba2 <- conf_fcast |>
          filter(grepl("BA\\.2|BA2", lineage, ignore.case = TRUE) |
                   lineage == ba2_lineage[1])

        if (nrow(conf_ba2) > 0) {
          conformal_trigger <- conf_ba2 |>
            filter(lower_95 > TRIGGER_THRESHOLD) |>
            slice_min(origin_date) |>
            pull(origin_date) |>
            min(na.rm = TRUE)
        }
      }
    }

    # Recalibrated trigger
    recal_trigger <- NA
    if (!is.null(conformal_cal$recalibrated)) {
      recal_fcast <- if (is.data.frame(conformal_cal$recalibrated)) {
        conformal_cal$recalibrated
      } else {
        conformal_cal$recalibrated$forecasts
      }

      if (!is.null(recal_fcast)) {
        recal_ba2 <- recal_fcast |>
          filter(grepl("BA\\.2|BA2", lineage, ignore.case = TRUE) |
                   lineage == ba2_lineage[1])

        if (nrow(recal_ba2) > 0) {
          recal_trigger <- recal_ba2 |>
            filter(lower_95 > TRIGGER_THRESHOLD) |>
            slice_min(origin_date) |>
            pull(origin_date) |>
            min(na.rm = TRUE)
        }
      }
    }

    # Compute trigger errors
    trigger_df <- tibble(
      dataset = bt$dataset,
      method  = c("Parametric", "Conformal", "Recalibrated"),
      trigger_date = c(parametric_trigger, conformal_trigger, recal_trigger),
      actual_crossing_date = actual_crossing,
      error_days = as.numeric(
        c(parametric_trigger, conformal_trigger, recal_trigger) - actual_crossing
      )
    )

    # False trigger rate: proportion of origins where method triggered
    # but actual hadn't crossed yet
    compute_false_rate <- function(fcast_df, lineage_name) {
      if (is.null(fcast_df) || nrow(fcast_df) == 0) return(NA_real_)
      lin_data <- fcast_df |>
        filter(lineage == lineage_name)
      if (nrow(lin_data) == 0) return(NA_real_)

      triggered <- lin_data |>
        filter(lower_95 > TRIGGER_THRESHOLD)
      false_triggers <- triggered |>
        filter(observed <= TRIGGER_THRESHOLD)

      if (nrow(triggered) == 0) return(0)
      nrow(false_triggers) / nrow(triggered)
    }

    trigger_df$false_trigger_rate <- c(
      compute_false_rate(ba2_fcast, ba2_lineage[1]),
      tryCatch(compute_false_rate(
        if (is.data.frame(conformal_cal$conformal)) conformal_cal$conformal
        else conformal_cal$conformal$forecasts,
        ba2_lineage[1]
      ), error = function(e) NA_real_),
      tryCatch(compute_false_rate(
        if (is.data.frame(conformal_cal$recalibrated)) conformal_cal$recalibrated
        else conformal_cal$recalibrated$forecasts,
        ba2_lineage[1]
      ), error = function(e) NA_real_)
    )

    trigger_results[[key]] <- trigger_df
    cat(sprintf("      Actual crossing: %s\n", actual_crossing))
    cat(sprintf("      Parametric trigger: %s (%+d days)\n",
                parametric_trigger, trigger_df$error_days[1]))
    cat(sprintf("      Conformal trigger:  %s (%+d days)\n",
                conformal_trigger, trigger_df$error_days[2]))
    cat(sprintf("      Recalibrated trigger: %s (%+d days)\n",
                recal_trigger, trigger_df$error_days[3]))
  },
  error = function(e) {
    cat(sprintf("    WARNING: Trigger analysis failed for %s: %s\n",
                key, e$message))
  })
}

trigger_all <- bind_rows(trigger_results)

saveRDS(trigger_all, "analysis/results/decision_impact.rds")
cat("    Saved decision_impact.rds\n")

###############################################################################
# Section B: Sample size vs coverage curve
###############################################################################

cat("  [B] Sample size vs coverage analysis...\n")

# Systematically downsample the built-in BA.2 data to assess how sample size
# affects both point forecast accuracy (MAE) and interval calibration (coverage).
# This addresses the practical question: how many sequences per week are needed?

sample_sizes <- c(5000L, 2000L, 1000L, 500L, 200L, 100L)
n_replicates <- 50L

# Load built-in BA.2 data
ba2_data <- tryCatch(cdc_sarscov2_ba2, error = function(e) {
  cat("    WARNING: Built-in BA.2 data not available; skipping\n")
  NULL
})

if (!is.null(ba2_data)) {
  # Get original counts matrix
  original_counts <- ba2_data$counts
  original_dates  <- ba2_data$dates

  # Total sequences per time point
  original_totals <- rowSums(original_counts)
  # Original proportions
  original_props  <- original_counts / original_totals

  cat(sprintf("    Original data: %d time points, %d-%d seqs/period\n",
              nrow(original_counts),
              min(original_totals), max(original_totals)))

  # Run downsampling in parallel
  downsample_configs <- expand.grid(
    sample_size = sample_sizes,
    replicate   = seq_len(n_replicates),
    stringsAsFactors = FALSE
  )

  cat(sprintf("    Running %d downsampling configs (%d sizes × %d reps)...\n",
              nrow(downsample_configs), length(sample_sizes), n_replicates))

  run_downsample <- function(ss, rep_id, orig_props, orig_dates) {
    set.seed(20240414 + rep_id * 1000 + ss)

    # Resample counts from multinomial with reduced total
    n_times    <- nrow(orig_props)
    n_lineages <- ncol(orig_props)
    new_counts <- matrix(0L, nrow = n_times, ncol = n_lineages)

    for (t in seq_len(n_times)) {
      probs <- orig_props[t, ]
      probs[probs < 0] <- 0
      probs <- probs / sum(probs)
      new_counts[t, ] <- rmultinom(1, size = ss, prob = probs)
    }

    colnames(new_counts) <- colnames(orig_props)

    # Create lfq_data and run backtest
    tryCatch({
      ds_data <- lfq_data(counts = new_counts, dates = orig_dates)
      bt <- backtest(ds_data, engine = "mlr", horizons = c(14L),
                     min_train = 42L)

      fcast <- bt$forecasts

      mae <- mean(abs(fcast$predicted - fcast$observed), na.rm = TRUE)
      cov95 <- mean(fcast$observed >= fcast$lower_95 &
                       fcast$observed <= fcast$upper_95, na.rm = TRUE)

      tibble(
        sample_size = ss,
        replicate   = rep_id,
        mae         = mae,
        coverage_95 = cov95
      )
    },
    error = function(e) {
      tibble(sample_size = ss, replicate = rep_id,
             mae = NA_real_, coverage_95 = NA_real_)
    })
  }

  sample_size_results <- future_pmap_dfr(
    list(
      ss     = downsample_configs$sample_size,
      rep_id = downsample_configs$replicate
    ),
    run_downsample,
    orig_props = original_props,
    orig_dates = original_dates,
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )

  # Summarize by sample size
  sample_size_summary <- sample_size_results |>
    group_by(sample_size) |>
    summarise(
      mae_mean     = mean(mae, na.rm = TRUE),
      mae_sd       = sd(mae, na.rm = TRUE),
      mae_median   = median(mae, na.rm = TRUE),
      cov95_mean   = mean(coverage_95, na.rm = TRUE),
      cov95_sd     = sd(coverage_95, na.rm = TRUE),
      meets_mae5   = mean(mae < 0.05, na.rm = TRUE),      # MAE < 5%
      meets_cov90  = mean(coverage_95 > 0.90, na.rm = TRUE), # Coverage > 90%
      n_success    = sum(!is.na(mae)),
      .groups = "drop"
    ) |>
    arrange(desc(sample_size))

  cat("    Sample size analysis summary:\n")
  print(sample_size_summary |>
          select(sample_size, mae_mean, cov95_mean, meets_mae5, meets_cov90))

  saveRDS(list(raw = sample_size_results, summary = sample_size_summary),
          "analysis/results/sample_size_analysis.rds")
  cat("    Saved sample_size_analysis.rds\n")
} else {
  cat("    Skipped sample size analysis (no data)\n")
}

###############################################################################
# Section C: Training window length vs PIT uniformity
###############################################################################

cat("  [C] Training window vs PIT uniformity...\n")

# Hypothesis: Shorter training windows produce more uniform PIT distributions
# because the model assumption of constant growth advantage δ is less violated
# over shorter periods. This empirically supports the theoretical argument that
# model misspecification (time-varying δ) is the primary cause of miscalibration.

window_lengths <- c(4L, 6L, 8L, 10L, 12L, 14L)  # weeks
window_days    <- window_lengths * 7L

if (!is.null(ba2_data)) {
  window_results <- list()

  for (wl in window_days) {
    cat(sprintf("    Training window = %d days (%d weeks)...\n", wl, wl / 7))

    tryCatch({
      bt <- backtest(ba2_data, engine = "mlr",
                     horizons = 14L, min_train = wl)

      cal <- calibrate(bt)
      pit <- cal$pit_values
      ks  <- ks.test(pit, "punif")

      window_results[[as.character(wl)]] <- tibble(
        window_weeks = wl / 7,
        window_days  = wl,
        ks_D         = ks$statistic,
        ks_p         = ks$p.value,
        n_pit        = length(pit),
        pit_mean     = mean(pit),    # should be ~0.5 if uniform
        pit_sd       = sd(pit)       # should be ~0.289 if uniform (1/sqrt(12))
      )

      cat(sprintf("      KS D = %.3f, p = %.3g\n", ks$statistic, ks$p.value))
    },
    error = function(e) {
      cat(sprintf("      WARNING: Failed for window %d days: %s\n", wl, e$message))
    })
  }

  window_df <- bind_rows(window_results)

  # Expected pattern: KS D should decrease (better uniformity) with shorter windows
  if (nrow(window_df) >= 2) {
    trend <- cor(window_df$window_weeks, window_df$ks_D, method = "spearman")
    cat(sprintf("    Window-KS correlation (Spearman): %.3f\n", trend))
    if (trend > 0) {
      cat("    CONFIRMED: longer windows → worse calibration (higher KS D)\n")
      cat("    Supports: model misspecification accumulates with window length\n")
    } else {
      cat("    NOTE: Unexpected negative correlation; reporting honestly\n")
    }
  }

  saveRDS(window_df, "analysis/results/window_analysis.rds")
  cat("    Saved window_analysis.rds\n")
} else {
  cat("    Skipped window analysis (no data)\n")
}

cat("[04_decision_impact] Complete.\n")
