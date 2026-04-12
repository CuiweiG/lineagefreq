#!/usr/bin/env Rscript
# ============================================================
# run_all.R — Master script: run all analysis in order
# ============================================================
#
# Usage:
#   setwd("path/to/lineagefreq")
#   source("analysis/run_all.R")
#
# Expected runtime on AMD EPYC 9654 (192 cores): ~2 minutes
# Expected runtime on standard workstation (8 cores): ~10 minutes
# Expected runtime on laptop (4 cores): ~20 minutes
#
# Each script is independent after 00_setup.R has been run.
# Intermediate results are saved to analysis/results/, so
# individual scripts can be re-run without recomputing everything.
# ============================================================

total_start <- Sys.time()

cat("============================================================\n")
cat("lineagefreq v0.5.1 — Full validation analysis\n")
cat(sprintf("Started: %s\n", total_start))
cat("============================================================\n\n")

timings <- list()

run_section <- function(script_name) {
  t0 <- Sys.time()
  source(file.path("analysis", script_name), local = FALSE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timings[[script_name]] <<- elapsed
  cat(sprintf("  [%s: %.1f seconds]\n\n", script_name, elapsed))
}

run_section("00_setup.R")
run_section("01_benchmark.R")
run_section("02_calibration.R")
run_section("03_surveillance.R")
run_section("04_fitness.R")
run_section("05_influenza.R")
run_section("06_figures.R")
run_section("07_tables.R")

# ---- Summary ----
total_elapsed <- as.numeric(difftime(Sys.time(), total_start,
                                      units = "secs"))

cat("============================================================\n")
cat("Analysis complete.\n")
cat(sprintf("Total runtime: %.1f seconds (%.1f minutes)\n",
            total_elapsed, total_elapsed / 60))
cat("\nPer-section timing:\n")
for (nm in names(timings)) {
  cat(sprintf("  %-20s %6.1f s\n", nm, timings[[nm]]))
}
cat("\nOutput manifest:\n")
cat(sprintf("  Figures: %d PDF + %d PNG\n",
  length(list.files("analysis/figures", pattern = "\\.pdf$")),
  length(list.files("analysis/figures", pattern = "\\.png$"))))
cat(sprintf("  Tables:  %d .tex files\n",
  length(list.files("analysis/tables", pattern = "\\.tex$"))))
cat(sprintf("  Results: %d .rds files\n",
  length(list.files("analysis/results", pattern = "\\.rds$"))))
cat("============================================================\n")

# Clean up parallel workers
plan(sequential)
