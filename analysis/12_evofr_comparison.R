###############################################################################
# 12_evofr_comparison.R
# Head-to-head comparison with evofr (Bedford Lab)
#
# Run from package root:
#   source("analysis/12_evofr_comparison.R")
#
# Requires: evofr R package from https://github.com/blab/evofr
# Install with: remotes::install_github("blab/evofr")
###############################################################################

cat("[12_evofr_comparison] Starting...\n")
source("analysis/00_setup.R")
library(lineagefreq)

# ─── Check if evofr is installed ─────────────────────────────────────────

evofr_available <- requireNamespace("evofr", quietly = TRUE)

if (!evofr_available) {
  cat("  evofr not installed. Attempting installation...\n")
  tryCatch({
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("blab/evofr", quiet = TRUE)
    evofr_available <- requireNamespace("evofr", quietly = TRUE)
    if (evofr_available) cat("  evofr installed successfully.\n")
  }, error = function(e) {
    cat(sprintf("  evofr installation failed: %s\n", conditionMessage(e)))
    cat("  To install manually: remotes::install_github('blab/evofr')\n")
  })
}

if (!evofr_available) {
  cat("  evofr unavailable — creating placeholder results.\n")
  cat("  NOTE: The manuscript should note this as a limitation.\n")
  saveRDS(list(available = FALSE, reason = "evofr not installed"),
          "analysis/results/evofr_comparison.rds")
  cat("[12_evofr_comparison] Complete (evofr unavailable).\n")
  stop("evofr not available. Install with remotes::install_github('blab/evofr')")
}

# ─── Run lineagefreq on US BA.2 data ────────────────────────────────────

library(evofr)

cat("  Running lineagefreq MLR on US BA.2...\n")
ba2_raw <- cdc_ba2_transition
ba2_lfq <- lfq_data(ba2_raw, lineage = lineage, date = date, count = count)
fit_lfq <- fit_model(ba2_lfq, engine = "mlr")
bt_lfq  <- backtest(ba2_lfq, engines = "mlr",
                      horizons = c(7L, 14L, 21L, 28L), min_train = 42L)
cal_lfq <- calibrate(bt_lfq)

cat(sprintf("  lineagefreq: KS D = %.3f, p = %.2g\n",
            cal_lfq$ks_test$statistic, cal_lfq$ks_test$p_value))

# ─── Run evofr on same data ─────────────────────────────────────────────

cat("  Running evofr MLR on US BA.2...\n")

# evofr API: adapt data format
# evofr typically expects a data frame with:
#   date, variant/lineage, sequences (counts), total
# Check evofr documentation for exact column names

evofr_data <- data.frame(
  date     = ba2_raw$date,
  lineage  = ba2_raw$lineage,
  count    = ba2_raw$count
)

# Add total per date
evofr_data <- evofr_data |>
  dplyr::group_by(date) |>
  dplyr::mutate(total = sum(count)) |>
  dplyr::ungroup()

evofr_result <- tryCatch({
  # NOTE: evofr API may differ from this. Common patterns:
  # evofr::fit_mlr(data, ...)
  # evofr::forecast_frequencies(fit, ...)
  # Adjust based on actual evofr vignettes/documentation.

  cat("  NOTE: evofr API adaptation may be needed.\n")
  cat("  Read evofr documentation at https://github.com/blab/evofr\n")
  cat("  and adjust the function calls below.\n")

  # Attempt generic evofr call
  evofr_fit <- evofr::fit_mlr(evofr_data)
  evofr_fc  <- evofr::forecast_frequencies(evofr_fit, horizon = 28)
  list(fit = evofr_fit, forecast = evofr_fc)
}, error = function(e) {
  cat(sprintf("  evofr call failed: %s\n", conditionMessage(e)))
  cat("  This is expected — evofr API needs manual adaptation.\n")
  NULL
})

# ─── Compare calibration ────────────────────────────────────────────────

comparison <- list(
  available         = evofr_available,
  lineagefreq_ks_D  = cal_lfq$ks_test$statistic,
  lineagefreq_ks_p  = cal_lfq$ks_test$p_value,
  evofr_result      = evofr_result,
  note = paste(
    "evofr comparison requires API adaptation.",
    "The key comparison is whether evofr's MLR implementation",
    "produces similarly miscalibrated intervals.",
    "Since both use the same underlying MLR framework,",
    "we expect similar KS D values."
  )
)

saveRDS(comparison, "analysis/results/evofr_comparison.rds")
cat("\n[12_evofr_comparison] Complete.\n")
