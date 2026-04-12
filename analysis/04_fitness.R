#!/usr/bin/env Rscript
# ============================================================
# 04_fitness.R — Fitness decomposition analysis
# ============================================================
# Constructs an immune landscape for the US BA.2 transition
# period using published vaccination coverage data, then
# decomposes variant growth advantages into intrinsic
# transmissibility and immune escape components.
#
# Reference: Lyngse et al. (2022) found that BA.2 had primarily
# a transmissibility advantage over BA.1, not immune escape.
# ============================================================

cat("[04_fitness] Starting...\n")

if (!exists("theme_nature")) source("analysis/00_setup.R")

# ---- Fit model on BA.2 data ----
data(cdc_ba2_transition)
x_ba2 <- lfq_data(cdc_ba2_transition, lineage = lineage,
                   date = date, count = count)
fit_ba2 <- fit_model(x_ba2, engine = "mlr")

cat("[04_fitness] BA.2 growth rates:\n")
print(growth_advantage(fit_ba2))

# ---- Construct immunity landscape ----
# US vaccination data from Our World in Data (approximate):
#   Dec 2021: 62% primary series, 20% boosted
#   Jan 2022: 63% primary series, 25% boosted
#   Mar 2022: 65% primary series, 29% boosted
#   Jun 2022: 67% primary series, 32% boosted
#
# Immunity against each variant depends on both vaccination and
# prior infection. By Dec 2021, a large wave of Delta had occurred.
# The Omicron BA.1 wave began in Dec 2021 and peaked in Jan 2022.
#
# Estimated neutralising immunity (rough approximation):
#   BA.1: High immunity by Feb 2022 (massive BA.1 wave + vaccines)
#   BA.2: Moderate escape from BA.1 immunity (~30% reduced)
#   BA.2.12.1: Similar to BA.2 but slightly more escape
#   BA.4/5: Substantial escape from both BA.1 and vaccine immunity
#   Other (pre-Omicron): Very high population immunity
#
# These are deliberately rough. The decomposition result should
# be treated as a structured hypothesis, not a precise measurement.

cat("[04_fitness] Constructing immune landscape...\n")

dates_ba2 <- sort(unique(cdc_ba2_transition$date))
lineages  <- unique(cdc_ba2_transition$lineage)

# Immunity values per lineage (time-constant for simplicity)
# Higher immunity = more of the population is protected against this variant
imm_map <- c(
  "BA.1"      = 0.55,  # Widespread BA.1 infection + vaccination
  "BA.2"      = 0.35,  # Partial escape from BA.1 immunity
  "BA.2.12.1" = 0.30,  # Greater escape than BA.2
  "BA.4/5"    = 0.20,  # Substantial escape
  "Other"     = 0.65   # Pre-Omicron: high cumulative immunity
)

# Expand to all date × lineage combinations
imm_data <- expand.grid(
  date    = dates_ba2,
  lineage = lineages,
  stringsAsFactors = FALSE
) |>
  mutate(immunity = imm_map[lineage])

il <- immune_landscape(imm_data, date = date,
                       lineage = lineage, immunity = immunity)

# ---- Run decomposition ----
# Generation time: 3.2 days for Omicron BA.* subvariants
# (Du et al. 2022, Emerging Infectious Diseases)
cat("[04_fitness] Running fitness decomposition...\n")
fd <- fitness_decomposition(fit_ba2, il, generation_time = 3.2)
print(fd)

# ---- Validate against literature ----
# Lyngse et al. 2022 (Nature Communications):
# BA.2 had ~34% higher household SAR than BA.1, attributable
# primarily to higher intrinsic transmissibility, not immune escape.
#
# Our decomposition should show BA.2 with a high transmissibility
# fraction. If it shows mostly immune escape, that would be
# inconsistent with the published evidence and we should note this.

d_ba2 <- fd$decomposition[fd$decomposition$lineage == "BA.2", ]
if (nrow(d_ba2) > 0 && !is.na(d_ba2$transmissibility_fraction)) {
  cat(sprintf("\n[04_fitness] BA.2 decomposition:\n"))
  cat(sprintf("  Intrinsic: %.1f%%  Escape: %.1f%%\n",
    d_ba2$transmissibility_fraction * 100,
    d_ba2$escape_fraction * 100))

  if (d_ba2$transmissibility_fraction > 0.5) {
    cat("  Consistent with Lyngse et al. 2022: BA.2 advantage\n")
    cat("  is primarily intrinsic transmissibility.\n")
  } else {
    cat("  NOTE: Result suggests immune escape dominates.\n")
    cat("  This may reflect inaccurate immunity estimates.\n")
    cat("  The decomposition is sensitive to immunity inputs.\n")
  }
} else {
  cat("[04_fitness] BA.2 is the pivot lineage; no decomposition available.\n")
  cat("  The pivot is chosen as the highest-count lineage at the\n")
  cat("  first time point, which may be BA.1 or Other.\n")
}

# ---- Save results ----
saveRDS(list(
  fd = fd,
  il = il,
  fit = fit_ba2
), "analysis/results/fitness_results.rds")

cat("[04_fitness] Complete.\n\n")
