###############################################################################
# 07_influenza.R — Multi-pathogen demonstration (influenza)
# lineagefreq validation analysis
###############################################################################

cat("[07_influenza] Starting...\n")
source("analysis/00_setup.R")

set.seed(20240417)

flunet_data <- readRDS("analysis/results/flunet_prepared.rds")
data_source <- attr(flunet_data, "source")

cat(sprintf("  FluNet data source: %s\n", data_source))

influenza_results <- list()

if (length(flunet_data) > 0 && data_source == "who_flunet") {
  # ── Real FluNet data pipeline ────────────────────────────────────────────────
  cat("  Processing real FluNet data...\n")

  for (nm in names(flunet_data)) {
    cat(sprintf("    %s...\n", nm))

    flu_obj <- flunet_data[[nm]]

    tryCatch({
      flu_fit <- fit_model(flu_obj, engine = "mlr")

      flu_fcast <- tryCatch({
        forecast(flu_fit, horizon = 14L)
      }, error = function(e) {
        cat(sprintf("      WARNING: forecast() failed: %s\n", e$message))
        NULL
      })

      # backtest() param is 'engines' (plural)
      flu_bt <- tryCatch({
        backtest(flu_obj, engines = "mlr", horizons = c(7L, 14L),
                 min_train = 28L)
      }, error = function(e) {
        cat(sprintf("      WARNING: backtest() failed: %s\n", e$message))
        NULL
      })

      influenza_results[[nm]] <- list(
        data     = flu_obj,
        fit      = flu_fit,
        forecast = flu_fcast,
        backtest = flu_bt,
        source   = "who_flunet"
      )

      if (!is.null(flu_bt)) {
        # flu_bt is an lfq_backtest tibble; columns: predicted, observed, etc.
        mae <- mean(abs(flu_bt$predicted - flu_bt$observed), na.rm = TRUE)
        cat(sprintf("      MAE (14d): %.3f\n", mae))
      }
    },
    error = function(e) {
      cat(sprintf("      WARNING: Pipeline failed for %s: %s\n", nm, e$message))
    })
  }
} else {
  # ── Simulated data fallback ──────────────────────────────────────────────────
  cat("  Using simulated influenza data (real FluNet data unavailable).\n")
  cat("  NOTE: Results labeled as SIMULATED — not from real surveillance.\n")

  simulate_flu_season <- function(n_weeks = 52, total_per_week = 500,
                                  season_label = "2023") {
    weeks <- seq_len(n_weeks)
    dates <- as.Date(paste0(season_label, "-01-01")) + (weeks - 1) * 7

    intensity <- dnorm(weeks, mean = 10, sd = 6)
    intensity <- intensity / max(intensity)

    t_norm <- weeks / n_weeks
    p_h3   <- 0.5 * (1 - t_norm) + 0.15
    p_h1   <- 0.2 * t_norm + 0.10
    p_bvic <- 0.15 + 0.1 * sin(2 * pi * t_norm)
    p_byam <- 0.05 + 0.05 * cos(2 * pi * t_norm)

    p_mat <- cbind(p_h1, p_h3, p_bvic, p_byam)
    p_mat <- p_mat / rowSums(p_mat)
    p_mat[p_mat < 0] <- 0.01
    p_mat <- p_mat / rowSums(p_mat)

    n_per_week <- round(total_per_week * intensity)
    n_per_week <- pmax(n_per_week, 10)

    lineage_names <- c("AH1N1pdm09", "AH3N2", "B_Victoria", "B_Yamagata")

    # Build long-format data frame for lfq_data()
    season_weeks <- 1:26
    rows <- list()
    for (t in season_weeks) {
      probs <- p_mat[t, ]
      counts <- rmultinom(1, size = n_per_week[t], prob = probs)[, 1]
      for (j in seq_along(lineage_names)) {
        rows[[length(rows) + 1]] <- tibble(
          date    = dates[t],
          lineage = lineage_names[j],
          count   = counts[j]
        )
      }
    }
    sim_df <- bind_rows(rows)

    tryCatch({
      lfq_data(sim_df, lineage = lineage, date = date, count = count)
    },
    error = function(e) {
      cat(sprintf("    WARNING: lfq_data failed for simulated data: %s\n",
                  e$message))
      sim_df  # return raw data frame as fallback
    })
  }

  flu_sim_2023 <- simulate_flu_season(season_label = "2023")
  flu_sim_2024 <- simulate_flu_season(season_label = "2024")

  for (season_nm in c("flu_sim_2023", "flu_sim_2024")) {
    flu_obj <- get(season_nm)
    cat(sprintf("    %s...\n", season_nm))

    tryCatch({
      flu_fit <- fit_model(flu_obj, engine = "mlr")

      flu_fcast <- tryCatch(
        forecast(flu_fit, horizon = 14L),
        error = function(e) NULL
      )

      # backtest() param is 'engines' (plural)
      flu_bt <- tryCatch(
        backtest(flu_obj, engines = "mlr", horizons = c(7L, 14L),
                 min_train = 28L),
        error = function(e) NULL
      )

      influenza_results[[season_nm]] <- list(
        data     = flu_obj,
        fit      = flu_fit,
        forecast = flu_fcast,
        backtest = flu_bt,
        source   = "simulated"
      )

      if (!is.null(flu_bt)) {
        # flu_bt is the lfq_backtest tibble directly
        mae <- mean(abs(flu_bt$predicted - flu_bt$observed), na.rm = TRUE)
        cat(sprintf("      MAE (14d): %.3f\n", mae))
      }
    },
    error = function(e) {
      cat(sprintf("      WARNING: Pipeline failed: %s\n", e$message))
    })
  }
}

saveRDS(influenza_results, "analysis/results/influenza_results.rds")
cat(sprintf("  Saved influenza_results.rds (%d datasets)\n",
            length(influenza_results)))

cat("[07_influenza] Complete.\n")
