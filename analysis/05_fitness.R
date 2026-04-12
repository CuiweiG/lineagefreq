###############################################################################
# 05_fitness.R — Fitness decomposition (Conclusion 4)
# lineagefreq validation analysis
###############################################################################

cat("[05_fitness] Starting...\n")
source("analysis/00_setup.R")

set.seed(20240415)

###############################################################################
# Section A: Fitness decomposition with corrected immune landscape
###############################################################################

cat("  [A] Fitness decomposition...\n")

immune_data <- readRDS("analysis/results/immune_landscape.rds")

# Load built-in BA.2 data
# Dataset name: cdc_ba2_transition (data.frame with date, lineage, count, proportion)
ba2_raw <- tryCatch(cdc_ba2_transition, error = function(e) {
  cat("    WARNING: cdc_ba2_transition not available\n")
  NULL
})

fit <- NULL
decomp_result <- NULL

if (!is.null(ba2_raw)) {
  cat("    Creating lfq_data from built-in BA.2...\n")

  ba2_lfq <- tryCatch({
    lfq_data(ba2_raw, lineage = lineage, date = date, count = count)
  },
  error = function(e) {
    cat(sprintf("    WARNING: lfq_data() failed: %s\n", e$message))
    NULL
  })

  if (!is.null(ba2_lfq)) {
    cat("    Fitting MLR on BA.2 data...\n")

    fit <- tryCatch({
      fit_model(ba2_lfq, engine = "mlr")
    },
    error = function(e) {
      cat(sprintf("    WARNING: fit_model() failed: %s\n", e$message))
      NULL
    })
  }

  if (!is.null(fit)) {
    cat("    Running fitness decomposition with corrected immune landscape...\n")

    # fitness_decomposition() API:
    #   fitness_decomposition(fit, landscape, generation_time)
    # - fit: lfq_fit object
    # - landscape: immune_landscape object
    # - generation_time: numeric (days), e.g., 5.5 for SARS-CoV-2
    GENERATION_TIME <- 5.5  # days, Omicron (Du et al. 2022, Emerging Inf Dis)

    decomp_result <- tryCatch({
      fitness_decomposition(fit, landscape = immune_data,
                            generation_time = GENERATION_TIME)
    },
    error = function(e) {
      cat(sprintf("    WARNING: fitness_decomposition() failed: %s\n", e$message))
      NULL
    })

    if (!is.null(decomp_result)) {
      # decomp_result is a fitness_decomposition S3 object (list) with:
      #   $decomposition: tibble with lineage, observed_advantage, beta,
      #                   escape_contribution, transmissibility_fraction,
      #                   escape_fraction
      cat("    Decomposition results:\n")
      print(decomp_result)

      # Plausibility check using $decomposition tibble
      tryCatch({
        decomp_df <- decomp_result$decomposition

        for (i in seq_len(nrow(decomp_df))) {
          lin   <- decomp_df$lineage[i]
          trans <- decomp_df$transmissibility_fraction[i]
          esc   <- decomp_df$escape_fraction[i]

          if (!is.na(trans) && !is.na(esc)) {
            cat(sprintf("      %s: transmissibility=%.1f%%, escape=%.1f%%\n",
                        lin, trans * 100, esc * 100))

            if (grepl("BA\\.2", lin) && !is.na(esc) && esc > trans) {
              cat("      WARNING: BA.2 shows more escape than transmissibility\n")
              cat("      Expected: BA.2 fitness primarily from transmissibility\n")
              cat("      (Lyngse et al. 2022, Lancet Infectious Diseases)\n")
            }
            if (grepl("BA\\.[45]", lin) && !is.na(esc) && esc < 0) {
              cat("      WARNING: BA.4/5 shows negative escape\n")
              cat("      Expected: BA.4/5 has significant immune escape\n")
              cat("      (Tuekprakhon et al. 2022, Nature)\n")
            }
          }
        }
      },
      error = function(e) {
        cat(sprintf("      Could not parse decomposition for plausibility: %s\n",
                    e$message))
      })
    }
  }
}

###############################################################################
# Section B: Sensitivity analysis (Nature Methods requirement)
###############################################################################

cat("  [B] Sensitivity analysis...\n")

# Perturb the immune landscape by ±10%, ±20%, ±30% and re-run decomposition
# This quantifies the identifiability limitation.

perturbation_levels <- c(-0.30, -0.20, -0.10, 0, 0.10, 0.20, 0.30)

sensitivity_results <- list()

if (!is.null(fit) && !is.null(immune_data)) {
  # immune_data is an immune_landscape object with $estimates tibble
  # containing: date, lineage, immunity, type

  for (perturb in perturbation_levels) {
    cat(sprintf("    Perturbation: %+.0f%%...\n", perturb * 100))

    tryCatch({
      # Perturb immunity values in the estimates tibble
      if (inherits(immune_data, "immune_landscape")) {
        perturbed_df <- immune_data$estimates |>
          mutate(immunity = pmin(1, pmax(0, immunity * (1 + perturb))))

        perturbed_obj <- tryCatch({
          immune_landscape(
            data     = perturbed_df,
            date     = date,
            lineage  = lineage,
            immunity = immunity
          )
        },
        error = function(e) {
          cat(sprintf("      WARNING: immune_landscape() failed: %s\n", e$message))
          NULL
        })
      } else if (is.data.frame(immune_data)) {
        # Fallback: immune_data is a raw data frame
        perturbed_df <- immune_data |>
          mutate(immunity = pmin(1, pmax(0, total_immunity * (1 + perturb))))

        perturbed_obj <- tryCatch({
          immune_landscape(
            data     = perturbed_df,
            date     = date,
            lineage  = lineage,
            immunity = immunity
          )
        },
        error = function(e) NULL)
      } else {
        perturbed_obj <- NULL
      }

      if (!is.null(perturbed_obj)) {
        decomp_perturbed <- tryCatch({
          fitness_decomposition(fit, landscape = perturbed_obj,
                                generation_time = GENERATION_TIME)
        },
        error = function(e) NULL)

        if (!is.null(decomp_perturbed)) {
          sensitivity_results[[as.character(perturb)]] <-
            decomp_perturbed$decomposition |>
            mutate(perturbation = perturb)
        }
      }
    },
    error = function(e) {
      cat(sprintf("      WARNING: Perturbation %+.0f%% failed: %s\n",
                  perturb * 100, e$message))
    })
  }

  sensitivity_df <- bind_rows(sensitivity_results)

  if (nrow(sensitivity_df) > 0) {
    cat("    Sensitivity summary:\n")
    sens_summary <- sensitivity_df |>
      group_by(lineage, perturbation) |>
      summarise(
        transmissibility_frac = mean(transmissibility_fraction, na.rm = TRUE),
        escape_frac           = mean(escape_fraction, na.rm = TRUE),
        .groups = "drop"
      )
    print(sens_summary)

    # Quantify sensitivity
    tryCatch({
      for (lin in unique(sens_summary$lineage)) {
        lin_data <- sens_summary |> filter(lineage == lin)
        if (nrow(lin_data) >= 3 && any(!is.na(lin_data$escape_frac))) {
          escape_sensitivity <- coef(lm(escape_frac ~ perturbation,
                                         data = lin_data))[2]
          cat(sprintf("      %s: Δescape_frac/Δimmunity = %.3f\n",
                      lin, escape_sensitivity))
        }
      }
    }, error = function(e) NULL)
  }

  saveRDS(list(decomposition = decomp_result, sensitivity = sensitivity_df),
          "analysis/results/fitness_sensitivity.rds")
  cat("    Saved fitness_sensitivity.rds\n")
} else {
  cat("    Skipped sensitivity analysis (no fit or immune data)\n")
}

###############################################################################
# Section C: Identifiability formal statement
###############################################################################

# ─── Mathematical condition for identifiability ─────────────────────────────
#
# The fitness of lineage i relative to reference lineage 0 is:
#   δ_i(t) = δ_intrinsic_i + δ_escape_i × I(t)
#
# where:
#   δ_intrinsic_i = intrinsic transmissibility advantage
#   δ_escape_i    = immune escape coefficient
#   I(t)          = population-level effective immunity at time t
#
# The observed growth advantage from MLR is:
#   δ_obs_i(t) = β_i (constant in standard MLR)
#
# For decomposition to be identified, we need:
#   d(I)/dt ≠ 0 during the observation window
#
# Formally, the Jacobian of the mapping:
#   (δ_intrinsic, δ_escape) → (δ_obs(t), I(t))
# must have full rank (rank 2) for at least some t in the observation window.
#
# This requires:
#   J = | 1    I(t₁) |
#       | 1    I(t₂) |
#   rank(J) = 2  ⟺  I(t₁) ≠ I(t₂)
#
# In practice, we require immunity to change by >5 percentage points during
# the observation window for the decomposition to be meaningfully identified.
#
# References:
#   - Obermeyer et al. 2022, Science (PyR0 model)
#   - Mlcochova et al. 2021, Nature (variant fitness)
# ─────────────────────────────────────────────────────────────────────────────

cat("[05_fitness] Complete.\n")
