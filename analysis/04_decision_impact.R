###############################################################################
# 04_decision_impact.R — Decision impact analysis (NEW — Nature Methods)
# lineagefreq validation analysis
###############################################################################

cat("[04_decision_impact] Starting...\n")
source("analysis/00_setup.R")

set.seed(20240414)

benchmark   <- readRDS("analysis/results/benchmark_multicountry.rds")
calibration <- readRDS("analysis/results/calibration_comparison.rds")

###############################################################################
# Section A: Vaccine update trigger analysis
###############################################################################

cat("  [A] Vaccine update trigger analysis...\n")

# Scenario: A public health agency triggers a vaccine composition update review
# when a new variant's 95% LOWER confidence bound exceeds 30% frequency.
# We compare trigger timing across three interval estimation methods:
# - Parametric: standard MLR-based intervals (potentially too narrow)
# - Conformal: split conformal intervals from 03_calibration.R
# - Recalibrated: scaled parametric intervals from 03_calibration.R

TRIGGER_THRESHOLD <- 0.30  # 30% frequency

# Find BA.2 backtest results (US built-in data preferred)
ba2_keys <- names(benchmark$backtest_results) |>
  grep("BA2|ba2|ba2_transition", x = _, value = TRUE)

if (length(ba2_keys) == 0) {
  cat("    WARNING: No BA.2 backtest results found; using first available\n")
  ba2_keys <- names(benchmark$backtest_results)[1]
}

trigger_results <- list()

for (key in ba2_keys) {
  bt <- benchmark$backtest_results[[key]]

  tryCatch({
    cat(sprintf("    Analyzing trigger timing for %s...\n", bt$dataset))

    # backtest() returns an lfq_backtest tibble directly (no $forecasts)
    # Columns: origin_date, target_date, horizon, engine, lineage, predicted,
    #          lower, upper, observed
    bt_tibble <- bt$result

    # Identify BA.2 lineage
    ba2_lineage <- bt_tibble |>
      filter(grepl("BA\\.2|BA2", lineage, ignore.case = TRUE)) |>
      pull(lineage) |>
      unique()

    if (length(ba2_lineage) == 0) {
      cat("      No BA.2 lineage found; using lineage with largest max frequency\n")
      ba2_lineage <- bt_tibble |>
        group_by(lineage) |>
        summarise(max_obs = max(observed, na.rm = TRUE), .groups = "drop") |>
        slice_max(max_obs, n = 1) |>
        pull(lineage)
    }

    ba2_bt <- bt_tibble |>
      filter(lineage == ba2_lineage[1])

    # Actual crossing date: first origin_date when BA.2 observed > 30%
    actual_crossing <- ba2_bt |>
      filter(observed > TRIGGER_THRESHOLD) |>
      slice_min(origin_date, n = 1) |>
      pull(origin_date) |>
      min(na.rm = TRUE)

    # Parametric trigger: first origin where lower > threshold
    # Column is 'lower' (not 'lower_95')
    parametric_trigger <- ba2_bt |>
      filter(!is.na(lower), lower > TRIGGER_THRESHOLD) |>
      slice_min(origin_date, n = 1) |>
      pull(origin_date)
    parametric_trigger <- if (length(parametric_trigger) > 0) min(parametric_trigger) else NA

    # Conformal trigger: use conformal intervals from 03_calibration.R
    conformal_trigger <- NA
    cal_comp <- calibration$calibration_comparison[[key]]
    if (!is.null(cal_comp) && !is.null(cal_comp$conformal)) {
      conf_data <- cal_comp$conformal
      conf_ba2 <- conf_data |>
        filter(lineage == ba2_lineage[1])

      if (nrow(conf_ba2) > 0 && "conf_lower" %in% names(conf_ba2)) {
        triggered <- conf_ba2 |>
          filter(!is.na(conf_lower), conf_lower > TRIGGER_THRESHOLD) |>
          slice_min(origin_date, n = 1) |>
          pull(origin_date)
        conformal_trigger <- if (length(triggered) > 0) min(triggered) else NA
      }
    }

    # Recalibrated trigger: use recalibrated intervals from 03_calibration.R
    recal_trigger <- NA
    if (!is.null(cal_comp) && !is.null(cal_comp$recalibrated)) {
      recal_data <- cal_comp$recalibrated
      recal_ba2 <- recal_data |>
        filter(lineage == ba2_lineage[1])

      if (nrow(recal_ba2) > 0 && "recal_lower" %in% names(recal_ba2)) {
        triggered <- recal_ba2 |>
          filter(!is.na(recal_lower), recal_lower > TRIGGER_THRESHOLD) |>
          slice_min(origin_date, n = 1) |>
          pull(origin_date)
        recal_trigger <- if (length(triggered) > 0) min(triggered) else NA
      }
    }

    # Compute trigger errors
    trigger_df <- tibble(
      dataset = bt$dataset,
      method  = c("Parametric", "Conformal", "Recalibrated"),
      trigger_date = as.Date(c(parametric_trigger, conformal_trigger, recal_trigger)),
      actual_crossing_date = actual_crossing,
      error_days = as.numeric(
        as.Date(c(parametric_trigger, conformal_trigger, recal_trigger)) - actual_crossing
      )
    )

    # False trigger rate: proportion of origins where method triggered
    # but actual hadn't crossed yet
    compute_false_rate <- function(bt_df, lower_col, lineage_name) {
      lin_data <- bt_df |> filter(lineage == lineage_name)
      if (nrow(lin_data) == 0 || !lower_col %in% names(lin_data)) return(NA_real_)

      triggered <- lin_data |> filter(.data[[lower_col]] > TRIGGER_THRESHOLD)
      if (nrow(triggered) == 0) return(0)
      false_triggers <- triggered |> filter(observed <= TRIGGER_THRESHOLD)
      nrow(false_triggers) / nrow(triggered)
    }

    trigger_df$false_trigger_rate <- c(
      compute_false_rate(ba2_bt, "lower", ba2_lineage[1]),
      if (!is.null(cal_comp$conformal))
        compute_false_rate(cal_comp$conformal, "conf_lower", ba2_lineage[1])
      else NA_real_,
      if (!is.null(cal_comp$recalibrated))
        compute_false_rate(cal_comp$recalibrated, "recal_lower", ba2_lineage[1])
      else NA_real_
    )

    trigger_results[[key]] <- trigger_df
    cat(sprintf("      Actual crossing: %s\n", actual_crossing))
    for (i in 1:3) {
      cat(sprintf("      %s trigger: %s (%s days)\n",
                  trigger_df$method[i],
                  ifelse(is.na(trigger_df$trigger_date[i]), "Not triggered",
                         as.character(trigger_df$trigger_date[i])),
                  ifelse(is.na(trigger_df$error_days[i]), "N/A",
                         sprintf("%+d", as.integer(trigger_df$error_days[i])))))
    }
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

sample_sizes <- c(5000L, 2000L, 1000L, 500L, 200L, 100L)
n_replicates <- 50L

# Load built-in BA.2 data
# Dataset name is cdc_ba2_transition (not cdc_sarscov2_ba2)
ba2_raw <- tryCatch(cdc_ba2_transition, error = function(e) {
  cat("    WARNING: Built-in BA.2 data not available; skipping\n")
  NULL
})

if (!is.null(ba2_raw)) {
  # ba2_raw is a data.frame with columns: date, lineage, count, proportion
  # Extract unique dates and lineages to compute proportions
  ba2_wide <- ba2_raw |>
    group_by(date, lineage) |>
    summarise(count = sum(count), .groups = "drop") |>
    group_by(date) |>
    mutate(total = sum(count), prop = count / total) |>
    ungroup()

  unique_dates    <- sort(unique(ba2_wide$date))
  unique_lineages <- sort(unique(ba2_wide$lineage))

  # Build proportions matrix for downsampling
  prop_matrix <- ba2_wide |>
    select(date, lineage, prop) |>
    pivot_wider(names_from = lineage, values_from = prop, values_fill = 0) |>
    arrange(date)
  original_dates <- prop_matrix$date
  original_props <- as.matrix(prop_matrix |> select(-date))

  cat(sprintf("    Original data: %d time points, %d lineages\n",
              length(original_dates), ncol(original_props)))

  # Create lfq_data from built-in for backtest validation
  ba2_lfq <- tryCatch(
    lfq_data(ba2_raw, lineage = lineage, date = date, count = count),
    error = function(e) {
      cat(sprintf("    WARNING: lfq_data() failed: %s\n", e$message))
      NULL
    }
  )

  # Run downsampling in parallel
  downsample_configs <- expand.grid(
    sample_size = sample_sizes,
    replicate   = seq_len(n_replicates),
    stringsAsFactors = FALSE
  )

  cat(sprintf("    Running %d downsampling configs (%d sizes × %d reps)...\n",
              nrow(downsample_configs), length(sample_sizes), n_replicates))

  run_downsample <- function(ss, rep_id, orig_props, orig_dates, lineage_names) {
    set.seed(20240414 + rep_id * 1000 + ss)

    n_times <- nrow(orig_props)

    # Resample counts from multinomial with reduced total
    count_rows <- list()
    for (t in seq_len(n_times)) {
      probs <- orig_props[t, ]
      probs[probs < 0] <- 0
      probs <- probs / sum(probs)
      sampled <- rmultinom(1, size = ss, prob = probs)[, 1]

      for (j in seq_along(lineage_names)) {
        count_rows[[length(count_rows) + 1]] <- tibble(
          date = orig_dates[t], lineage = lineage_names[j], count = sampled[j]
        )
      }
    }
    ds_df <- bind_rows(count_rows)

    # Create lfq_data and run backtest
    tryCatch({
      ds_data <- lfq_data(ds_df, lineage = lineage, date = date, count = count)
      # backtest() param is 'engines' (plural), not 'engine'
      bt <- backtest(ds_data, engines = "mlr", horizons = 14L, min_train = 42L)

      # bt is an lfq_backtest tibble; columns: predicted, observed, lower, upper
      mae <- mean(abs(bt$predicted - bt$observed), na.rm = TRUE)
      cov95 <- mean(bt$observed >= bt$lower &
                     bt$observed <= bt$upper, na.rm = TRUE)

      tibble(sample_size = ss, replicate = rep_id, mae = mae, coverage_95 = cov95)
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
    orig_props     = original_props,
    orig_dates     = original_dates,
    lineage_names  = colnames(original_props),
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
      meets_mae5   = mean(mae < 0.05, na.rm = TRUE),
      meets_cov90  = mean(coverage_95 > 0.90, na.rm = TRUE),
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
# over shorter periods.

window_lengths <- c(4L, 6L, 8L, 10L, 12L, 14L)  # weeks
window_days    <- window_lengths * 7L

# Use the lfq_data object created from built-in BA.2 data
if (exists("ba2_lfq") && !is.null(ba2_lfq)) {
  window_results <- list()

  for (wl in window_days) {
    cat(sprintf("    Training window = %d days (%d weeks)...\n", wl, wl / 7))

    tryCatch({
      # backtest() param is 'engines' (plural)
      bt <- backtest(ba2_lfq, engines = "mlr", horizons = 14L, min_train = wl)

      # calibrate() accepts lfq_backtest directly
      cal <- calibrate(bt)
      pit <- cal$pit_values
      ks  <- ks.test(pit, "punif")

      window_results[[as.character(wl)]] <- tibble(
        window_weeks = wl / 7,
        window_days  = wl,
        ks_D         = ks$statistic,
        ks_p         = ks$p.value,
        n_pit        = length(pit),
        pit_mean     = mean(pit),
        pit_sd       = sd(pit)
      )

      cat(sprintf("      KS D = %.3f, p = %.3g\n", ks$statistic, ks$p.value))
    },
    error = function(e) {
      cat(sprintf("      WARNING: Failed for window %d days: %s\n", wl, e$message))
    })
  }

  window_df <- bind_rows(window_results)

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
  cat("    Skipped window analysis (no BA.2 lfq_data)\n")
}

cat("[04_decision_impact] Complete.\n")
