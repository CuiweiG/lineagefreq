###############################################################################
# 02_benchmark.R — Multi-country forecasting benchmark (Conclusion 1)
# lineagefreq validation analysis
###############################################################################

cat("[02_benchmark] Starting...\n")
source("analysis/00_setup.R")

set.seed(20240412)

###############################################################################
# Section A: Multi-country backtest
###############################################################################

cat("  [A] Multi-country backtest...\n")

ecdc_data <- readRDS("analysis/results/ecdc_prepared.rds")
cat(sprintf("    Loaded %d ECDC datasets\n", length(ecdc_data)))

# Also include package built-in datasets
builtin_datasets <- list()
tryCatch({
  builtin_datasets[["US_BA2_builtin"]] <- cdc_sarscov2_ba2
  cat("    Loaded built-in cdc_sarscov2_ba2\n")
}, error = function(e) cat(sprintf("    WARNING: cdc_sarscov2_ba2 not available: %s\n", e$message)))

tryCatch({
  builtin_datasets[["US_JN1_builtin"]] <- cdc_sarscov2_jn1
  cat("    Loaded built-in cdc_sarscov2_jn1\n")
}, error = function(e) cat(sprintf("    WARNING: cdc_sarscov2_jn1 not available: %s\n", e$message)))

all_datasets <- c(ecdc_data, builtin_datasets)

# Engines to test
engines <- c("mlr", "piantham")

# Backtest parameters
horizons  <- c(7L, 14L, 21L, 28L)
min_train <- 42L  # 6 weeks minimum training window

# Run backtests in parallel across datasets
# For each dataset × engine combination
backtest_configs <- expand.grid(
  dataset = names(all_datasets),
  engine  = engines,
  stringsAsFactors = FALSE
)

cat(sprintf("    Running %d backtest configurations across %d workers...\n",
            nrow(backtest_configs), nbrOfWorkers()))

run_single_backtest <- function(dataset_name, engine_name, datasets_list) {
  data_obj <- datasets_list[[dataset_name]]

  # Time the backtest
  timing <- system.time({
    result <- tryCatch({
      backtest(data_obj, engine = engine_name,
               horizons = horizons, min_train = min_train)
    },
    error = function(e) {
      # Return a plain list with $error on failure; successful backtest()
      # returns an lfq_backtest tibble (no $error field).
      list(error = conditionMessage(e))
    })
  })

  list(
    dataset = dataset_name,
    engine  = engine_name,
    result  = result,
    elapsed = timing["elapsed"],
    cores   = nbrOfWorkers()
  )
}

# Parallelize across datasets using furrr
backtest_results <- future_pmap(
  list(
    dataset_name = backtest_configs$dataset,
    engine_name  = backtest_configs$engine
  ),
  run_single_backtest,
  datasets_list = all_datasets,
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

names(backtest_results) <- paste(backtest_configs$dataset,
                                  backtest_configs$engine, sep = "__")

# Separate successful from failed
# backtest() returns an lfq_backtest tibble on success (inherits "lfq_backtest")
# and a plain list with $error on failure
is_failed <- function(x) is.list(x$result) && !is.null(x$result$error)
successful <- Filter(function(x) !is_failed(x), backtest_results)
failed     <- Filter(is_failed, backtest_results)

cat(sprintf("    Completed: %d successful, %d failed\n",
            length(successful), length(failed)))

for (f in failed) {
  cat(sprintf("    FAILED: %s/%s — %s\n", f$dataset, f$engine, f$result$error))
}

# Diagnostic: print structure of first successful result
if (length(successful) > 0) {
  first_bt <- successful[[1]]$result
  cat(sprintf("    Backtest result class: %s\n",
              paste(class(first_bt), collapse = ", ")))
  cat(sprintf("    Backtest result columns: %s\n",
              paste(names(first_bt), collapse = ", ")))
  cat(sprintf("    Backtest result rows: %d\n", nrow(first_bt)))
}

###############################################################################
# Compute accuracy metrics per dataset-engine
###############################################################################

cat("  Computing accuracy metrics...\n")

# Use score_forecasts() if available; otherwise compute manually.
# score_forecasts() returns a tibble with columns: engine, horizon, metric, value
compute_metrics <- function(bt_result) {
  tryCatch({
    # backtest() returns an lfq_backtest tibble directly (no $forecasts slot)
    bt <- bt_result$result

    # Try score_forecasts() first — it computes MAE, CRPS, coverage, etc.
    scores <- tryCatch({
      score_forecasts(bt)
    }, error = function(e) NULL)

    if (!is.null(scores)) {
      # score_forecasts() returns long format: engine, horizon, metric, value
      # Pivot to wide for our metrics table
      scores_wide <- scores |>
        select(horizon, metric, value) |>
        pivot_wider(names_from = metric, values_from = value)

      # Add median AE manually (score_forecasts may not include it)
      median_ae_by_h <- bt |>
        group_by(horizon) |>
        summarise(median_ae = median(abs(predicted - observed), na.rm = TRUE),
                  n_origins = n_distinct(origin_date),
                  .groups = "drop")

      result <- scores_wide |>
        left_join(median_ae_by_h, by = "horizon")

      return(result)
    }

    # Fallback: compute metrics manually
    metrics_list <- lapply(horizons, function(h) {
      horizon_data <- bt |> filter(horizon == h)
      if (nrow(horizon_data) == 0) return(NULL)

      mae       <- mean(abs(horizon_data$predicted - horizon_data$observed), na.rm = TRUE)
      median_ae <- median(abs(horizon_data$predicted - horizon_data$observed), na.rm = TRUE)

      # Coverage at 95%: backtest returns 'lower' and 'upper' columns (95% CI)
      coverage_95 <- if (all(c("lower", "upper") %in% names(horizon_data))) {
        mean(horizon_data$observed >= horizon_data$lower &
             horizon_data$observed <= horizon_data$upper, na.rm = TRUE)
      } else {
        NA_real_
      }

      tibble(
        horizon     = h,
        mae         = mae,
        median_ae   = median_ae,
        coverage_95 = coverage_95,
        n_origins   = n_distinct(horizon_data$origin_date)
      )
    })

    bind_rows(metrics_list)
  },
  error = function(e) {
    cat(sprintf("    WARNING: Metrics computation failed: %s\n", e$message))
    tibble()
  })
}

# Block bootstrap for MAE confidence intervals
# Block size = 3 origins (to preserve temporal correlation)
# Reference: Künsch 1989, Annals of Statistics
compute_bootstrap_ci <- function(bt_result, n_boot = 1000, block_size = 3) {
  tryCatch({
    # backtest() returns the tibble directly
    bt <- bt_result$result

    ci_list <- lapply(horizons, function(h) {
      horizon_data <- bt |> filter(horizon == h)

      if (nrow(horizon_data) < block_size * 2) {
        return(tibble(horizon = h, mae_lower = NA, mae_upper = NA))
      }

      errors <- abs(horizon_data$predicted - horizon_data$observed)
      n <- length(errors)

      # Block bootstrap: tsboot from boot package
      mae_stat <- function(tsb) mean(tsb)
      boot_result <- tryCatch(
        tsboot(errors, statistic = mae_stat, R = n_boot,
               l = block_size, sim = "fixed"),
        error = function(e) NULL
      )

      if (!is.null(boot_result)) {
        ci <- quantile(boot_result$t, c(0.025, 0.975), na.rm = TRUE)
        tibble(horizon = h, mae_lower = ci[1], mae_upper = ci[2])
      } else {
        tibble(horizon = h, mae_lower = NA, mae_upper = NA)
      }
    })

    bind_rows(ci_list)
  },
  error = function(e) {
    tibble()
  })
}

all_metrics <- list()
all_ci      <- list()

for (nm in names(successful)) {
  bt <- successful[[nm]]
  cat(sprintf("    Metrics for %s/%s...\n", bt$dataset, bt$engine))

  m <- compute_metrics(bt)
  if (nrow(m) > 0) {
    m$dataset <- bt$dataset
    m$engine  <- bt$engine
    all_metrics[[nm]] <- m
  }

  ci <- compute_bootstrap_ci(bt)
  if (nrow(ci) > 0) {
    ci$dataset <- bt$dataset
    ci$engine  <- bt$engine
    all_ci[[nm]] <- ci
  }
}

metrics_df <- bind_rows(all_metrics)
ci_df      <- bind_rows(all_ci)

cat(sprintf("    Metrics computed for %d configurations\n", length(all_metrics)))

###############################################################################
# Section B: Bedford Lab comparison table
###############################################################################

cat("  [B] Bedford Lab comparison...\n")

# Hardcoded reference values from Abousamra et al. 2024, PLOS Computational
# Biology, Table 1. These are for US CDC Nowcast data, MLR engine.
# Source: https://doi.org/10.1371/journal.pcbi.1012443
bedford_ref <- tibble(
  source    = "Bedford Lab (Abousamra et al. 2024)",
  engine    = "mlr",
  dataset   = "US CDC Nowcast",
  horizon   = c(7, 14, 30),
  median_ae = c(0.006, NA, NA),       # ~0.6% at 7d
  mae       = c(0.030, 0.060, 0.090)  # ~3%, ~6%, ~8-10% mean AE
)

# Merge with our results for matching horizons
comparison_df <- metrics_df |>
  filter(engine == "mlr") |>
  select(dataset, engine, horizon, mae, median_ae) |>
  mutate(source = "lineagefreq") |>
  bind_rows(bedford_ref)

cat("    Bedford Lab comparison table created\n")

###############################################################################
# Section C: Runtime benchmark
###############################################################################

cat("  [C] Runtime benchmark...\n")

runtime_df <- tibble(
  dataset    = sapply(successful, `[[`, "dataset"),
  engine     = sapply(successful, `[[`, "engine"),
  elapsed_s  = sapply(successful, function(x) as.numeric(x$elapsed)),
  cores_used = sapply(successful, `[[`, "cores")
) |>
  mutate(
    # backtest() stores n_origins as an attribute, not a column
    n_origins = sapply(names(successful), function(nm) {
      tryCatch(
        attr(successful[[nm]]$result, "n_origins"),
        error = function(e) {
          # Fallback: count distinct origin_date
          tryCatch(
            n_distinct(successful[[nm]]$result$origin_date),
            error = function(e2) NA_integer_
          )
        }
      )
    })
  )

cat("    Runtime summary:\n")
print(runtime_df |> select(dataset, engine, elapsed_s, cores_used))

###############################################################################
# Save all results
###############################################################################

benchmark_results <- list(
  backtest_results = successful,
  failed_results   = failed,
  metrics          = metrics_df,
  bootstrap_ci     = ci_df,
  bedford_ref      = bedford_ref,
  comparison       = comparison_df,
  runtime          = runtime_df,
  config = list(
    horizons  = horizons,
    min_train = min_train,
    engines   = engines,
    n_boot    = 1000,
    block_size = 3
  )
)

saveRDS(benchmark_results, "analysis/results/benchmark_multicountry.rds")
cat("    Saved benchmark_multicountry.rds\n")

cat("[02_benchmark] Complete.\n")
