###############################################################################
# 10_caps_validation.R — Validate CAPS against existing methods
# lineagefreq validation analysis
#
# Run from package root after devtools::install() or devtools::load_all():
#   source("analysis/10_caps_validation.R")
###############################################################################

cat("[10_caps_validation] Starting...\n")
source("analysis/00_setup.R")

library(lineagefreq)

# ─── Load data ────────────────────────────────────────────────────────────

ecdc  <- readRDS("analysis/results/ecdc_prepared.rds")
bench <- readRDS("analysis/results/benchmark_multicountry.rds")

datasets_to_test <- c("Denmark_BA2", "France_BA2", "Germany_BA2",
                       "Netherlands_BA2", "Spain_BA2")

# Built-in US BA.2 data
ba2_raw <- cdc_ba2_transition
ba2_lfq <- lfq_data(ba2_raw, lineage = lineage, date = date, count = count)

results_all    <- list()
comparison_all <- list()

# ─── Test 1: CAPS on built-in US BA.2 data ────────────────────────────

cat("\n[CAPS] Testing on US BA.2 data...\n")
fit_us <- fit_model(ba2_lfq, engine = "mlr")

for (method in c("caps", "caps_static", "caps_marginal")) {
  cat(sprintf("  Method: %s\n", method))
  caps_result <- tryCatch(
    caps_forecast(fit_us, ba2_lfq,
                  horizons = c(7L, 14L, 21L, 28L),
                  method = method, alpha = 0.05),
    error = function(e) {
      cat(sprintf("    FAILED: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(caps_result)) {
    print(caps_result)
    results_all[[paste0("US_BA2__", method)]] <- caps_result
  }
}

# ─── Test 2: CAPS on ECDC datasets ────────────────────────────────────

for (ds_name in datasets_to_test) {
  cat(sprintf("\n[CAPS] Testing on %s...\n", ds_name))

  ds <- ecdc[[ds_name]]
  if (is.null(ds)) {
    cat("  Dataset not found, skipping\n")
    next
  }

  fit <- tryCatch(
    fit_model(ds, engine = "mlr"),
    error = function(e) {
      cat(sprintf("  fit_model failed: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(fit)) next

  caps_result <- tryCatch(
    caps_forecast(fit, ds,
                  horizons = c(7L, 14L, 21L, 28L),
                  method = "caps", alpha = 0.05),
    error = function(e) {
      cat(sprintf("  CAPS failed: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(caps_result)) {
    cat(sprintf("  R_hat: %s\n",
                paste(sprintf("%.2f", caps_result$R_hat), collapse = ", ")))
    results_all[[paste0(ds_name, "__caps")]] <- caps_result
  }
}

# ─── Test 3: Evaluate CAPS vs split conformal vs parametric ───────────

cat("\n[CAPS] Generating comparison on US BA.2 backtest...\n")

# Run a fresh backtest for evaluation
bt_us <- backtest(ba2_lfq, engines = "mlr",
                   horizons = c(7L, 14L, 21L, 28L),
                   min_train = 42L)

caps_us <- results_all[["US_BA2__caps"]]

if (!is.null(caps_us)) {
  comp_df <- evaluate_caps(caps_us, bt_us)

  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("CAPS vs Split Conformal vs Parametric (US BA.2)\n")
  cat(strrep("=", 60), "\n")
  print(comp_df, width = 120)

  cat("\n-- CAPS advantage --\n")
  cat(sprintf("CAPS coverage >= 95%%: %d/%d horizons\n",
              sum(comp_df$caps_cov >= 0.95), nrow(comp_df)))
  cat(sprintf("CAPS width < conformal width: %d/%d horizons\n",
              sum(comp_df$caps_width < comp_df$conf_width), nrow(comp_df)))
  cat(sprintf("CAPS Winkler < conformal Winkler: %d/%d horizons\n",
              sum(comp_df$caps_winkler < comp_df$conf_winkler), nrow(comp_df)))

  comparison_all[["US_BA2"]] <- comp_df
}

# ─── Test 4: Evaluate CAPS on ECDC datasets ──────────────────────────

for (ds_name in datasets_to_test) {
  ds <- ecdc[[ds_name]]
  caps_result <- results_all[[paste0(ds_name, "__caps")]]
  if (is.null(ds) || is.null(caps_result)) next

  bt_ds <- tryCatch(
    backtest(ds, engines = "mlr",
             horizons = c(7L, 14L, 21L, 28L),
             min_train = 42L),
    error = function(e) NULL
  )

  if (!is.null(bt_ds)) {
    comp <- evaluate_caps(caps_result, bt_ds)
    cat(sprintf("\n[CAPS] %s comparison:\n", ds_name))
    print(comp)
    comparison_all[[ds_name]] <- comp
  }
}

# ─── Save results ────────────────────────────────────────────────────

if (length(comparison_all) > 0) {
  all_comp <- dplyr::bind_rows(comparison_all, .id = "dataset")
  saveRDS(all_comp, "analysis/results/caps_comparison.rds")
  cat("\nSaved analysis/results/caps_comparison.rds\n")

  # Summary across all datasets
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("CAPS Summary Across All Datasets\n")
  cat(strrep("=", 60), "\n")
  cat(sprintf("Total dataset-horizon pairs: %d\n", nrow(all_comp)))
  cat(sprintf("CAPS coverage >= 95%%: %d/%d (%.0f%%)\n",
              sum(all_comp$caps_cov >= 0.95), nrow(all_comp),
              100 * mean(all_comp$caps_cov >= 0.95)))
  cat(sprintf("CAPS beats conformal on Winkler: %d/%d (%.0f%%)\n",
              sum(all_comp$caps_winkler < all_comp$conf_winkler),
              nrow(all_comp),
              100 * mean(all_comp$caps_winkler < all_comp$conf_winkler)))
  cat(sprintf("Mean CAPS coverage: %.1f%%\n",
              100 * mean(all_comp$caps_cov, na.rm = TRUE)))
  cat(sprintf("Mean conformal coverage: %.1f%%\n",
              100 * mean(all_comp$conf_cov, na.rm = TRUE)))
  cat(sprintf("Mean parametric coverage: %.1f%%\n",
              100 * mean(all_comp$para_cov, na.rm = TRUE)))
}

saveRDS(results_all, "analysis/results/caps_results.rds")

cat("\n[10_caps_validation] Complete.\n")
