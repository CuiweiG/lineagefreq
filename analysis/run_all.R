###############################################################################
# run_all.R — Execute all analysis scripts in sequence
# lineagefreq validation analysis
#
# Expected total runtime on AMD EPYC 9654 (96 cores, 192 GB RAM):
#   ~20-40 minutes depending on data sizes and backtest complexity
#
# Usage: source("analysis/run_all.R") in RStudio or Rscript analysis/run_all.R
###############################################################################

cat("═══════════════════════════════════════════════════════════════\n")
cat("  lineagefreq — Nature Methods validation analysis\n")
cat("  Starting full pipeline\n")
cat(sprintf("  Time: %s\n", Sys.time()))
cat("═══════════════════════════════════════════════════════════════\n\n")

total_start <- proc.time()

scripts <- c(
  "analysis/00_setup.R",
  "analysis/01_data_prep.R",
  "analysis/02_benchmark.R",
  "analysis/03_calibration.R",
  "analysis/04_decision_impact.R",
  "analysis/05_fitness.R",
  "analysis/06_surveillance.R",
  "analysis/07_influenza.R",
  "analysis/08_figures.R",
  "analysis/09_tables.R"
)

timings <- list()

for (script in scripts) {
  script_name <- basename(script)
  cat(sprintf("\n{'='*60}\n  Running: %s\n{'='*60}\n", script_name))

  t0 <- proc.time()

  status <- tryCatch({
    source(script, local = FALSE)
    "SUCCESS"
  },
  error = function(e) {
    cat(sprintf("\n  *** ERROR in %s ***\n  %s\n\n", script_name,
                conditionMessage(e)))
    paste("FAILED:", conditionMessage(e))
  })

  elapsed <- (proc.time() - t0)["elapsed"]
  timings[[script_name]] <- list(status = status, elapsed = elapsed)

  cat(sprintf("\n  %s: %s (%.1f sec)\n", script_name, status, elapsed))
}

total_elapsed <- (proc.time() - total_start)["elapsed"]

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  Pipeline summary\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (nm in names(timings)) {
  t <- timings[[nm]]
  status_icon <- ifelse(t$status == "SUCCESS", "[OK]", "[!!]")
  cat(sprintf("  %s %-25s %6.1f sec  %s\n",
              status_icon, nm, t$elapsed, t$status))
}

n_ok   <- sum(sapply(timings, function(t) t$status == "SUCCESS"))
n_fail <- length(timings) - n_ok

cat(sprintf("\n  Total: %d/%d succeeded, %.1f sec (%.1f min)\n",
            n_ok, length(timings), total_elapsed, total_elapsed / 60))
cat(sprintf("  Completed: %s\n", Sys.time()))
cat("═══════════════════════════════════════════════════════════════\n")

if (n_fail > 0) {
  cat(sprintf("\n  WARNING: %d script(s) failed. Check output above.\n", n_fail))
}
