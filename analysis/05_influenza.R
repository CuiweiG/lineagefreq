#!/usr/bin/env Rscript
# ============================================================
# 05_influenza.R — Multi-pathogen demonstration
# ============================================================
# Demonstrates that lineagefreq works on non-SARS-CoV-2 data.
# Uses the included simulated influenza H3N2 dataset and also
# generates a custom seasonal influenza scenario.
# ============================================================

cat("[05_influenza] Starting...\n")

if (!exists("theme_nature")) source("analysis/00_setup.R")

# ---- Included H3N2 dataset ----
cat("[05_influenza] Analyzing included H3N2 dataset...\n")
data(influenza_h3n2)
x_flu <- lfq_data(influenza_h3n2, lineage = clade,
                  date = date, count = count)
fit_flu <- fit_model(x_flu, engine = "mlr")
fc_flu  <- forecast(fit_flu, horizon = 28)

cat("  Growth advantages (H3N2 clades):\n")
ga_flu <- growth_advantage(fit_flu)
print(ga_flu)

# ---- Custom seasonal scenario ----
# Simulate a realistic influenza season with 4 subtypes,
# one of which is replacing the others.
cat("[05_influenza] Simulating custom influenza scenario...\n")
sim_flu <- simulate_dynamics(
  n_lineages   = 4L,
  n_timepoints = 26L,   # 6-month season
  total_per_tp = 400L,   # typical sentinel surveillance volume
  advantages   = c(
    "H3N2_3C.2a1b" = 1.15,  # dominant, moderate growth
    "H1N1pdm09"    = 0.85,  # declining
    "B_Victoria"   = 0.70   # declining faster
  ),
  reference_name = "H3N2_3C.2a",
  start_date = as.Date("2024-10-01"),
  interval   = 7L,
  seed       = 2024
)

fit_sim_flu <- fit_model(sim_flu, engine = "mlr")
fc_sim_flu  <- forecast(fit_sim_flu, horizon = 28)

# ---- Save results ----
saveRDS(list(
  x_flu = x_flu,
  fit_flu = fit_flu,
  fc_flu = fc_flu,
  ga_flu = ga_flu,
  sim_flu = sim_flu,
  fit_sim_flu = fit_sim_flu,
  fc_sim_flu = fc_sim_flu
), "analysis/results/influenza_results.rds")

cat("[05_influenza] Complete.\n\n")
