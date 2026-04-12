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
      # Fit MLR model
      flu_fit <- fit_model(flu_obj, engine = "mlr")

      # Generate forecast
      flu_fcast <- tryCatch({
        forecast(flu_fit, horizon = 14)
      }, error = function(e) {
        cat(sprintf("      WARNING: forecast() failed: %s\n", e$message))
        NULL
      })

      # Run backtest
      flu_bt <- tryCatch({
        backtest(flu_obj, engine = "mlr", horizons = c(7L, 14L),
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
        mae <- mean(abs(flu_bt$forecasts$predicted - flu_bt$forecasts$observed),
                     na.rm = TRUE)
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

  # Simulate realistic influenza subtype dynamics
  # 4 subtypes: A(H1N1)pdm09, A(H3N2), B/Victoria, B/Yamagata
  # Simulate 2 seasons with typical patterns

  simulate_flu_season <- function(n_weeks = 52, total_per_week = 500,
                                  season_label = "2023") {
    weeks <- seq_len(n_weeks)
    dates <- as.Date(paste0(season_label, "-01-01")) + (weeks - 1) * 7

    # Seasonal intensity: peak around week 8-12 (Feb-Mar)
    intensity <- dnorm(weeks, mean = 10, sd = 6)
    intensity <- intensity / max(intensity)

    # Subtype proportions evolve over season
    # H3N2 dominant early, H1N1 rises later, B lineages minor
    t_norm <- weeks / n_weeks
    p_h3   <- 0.5 * (1 - t_norm) + 0.15
    p_h1   <- 0.2 * t_norm + 0.10
    p_bvic <- 0.15 + 0.1 * sin(2 * pi * t_norm)
    p_byam <- 0.05 + 0.05 * cos(2 * pi * t_norm)

    # Normalize
    p_mat <- cbind(p_h1, p_h3, p_bvic, p_byam)
    p_mat <- p_mat / rowSums(p_mat)
    p_mat[p_mat < 0] <- 0.01
    p_mat <- p_mat / rowSums(p_mat)

    # Generate counts
    n_per_week <- round(total_per_week * intensity)
    n_per_week <- pmax(n_per_week, 10)

    counts <- matrix(0L, nrow = n_weeks, ncol = 4)
    colnames(counts) <- c("AH1N1pdm09", "AH3N2", "B_Victoria", "B_Yamagata")

    for (t in seq_len(n_weeks)) {
      counts[t, ] <- rmultinom(1, size = n_per_week[t], prob = p_mat[t, ])
    }

    # Filter to flu season (weeks 1-26, roughly Oct-Mar)
    season_weeks <- 1:26
    counts <- counts[season_weeks, ]
    dates  <- dates[season_weeks]

    tryCatch({
      lfq_data(counts = counts, dates = dates)
    },
    error = function(e) {
      cat(sprintf("    WARNING: lfq_data failed for simulated data: %s\n",
                  e$message))
      list(counts = counts, dates = dates)
    })
  }

  # Generate two simulated seasons
  flu_sim_2023 <- simulate_flu_season(season_label = "2023")
  flu_sim_2024 <- simulate_flu_season(season_label = "2024")

  for (season_nm in c("flu_sim_2023", "flu_sim_2024")) {
    flu_obj <- get(season_nm)
    cat(sprintf("    %s...\n", season_nm))

    tryCatch({
      flu_fit <- fit_model(flu_obj, engine = "mlr")

      flu_fcast <- tryCatch(
        forecast(flu_fit, horizon = 14),
        error = function(e) NULL
      )

      flu_bt <- tryCatch(
        backtest(flu_obj, engine = "mlr", horizons = c(7L, 14L),
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
        mae <- mean(abs(flu_bt$forecasts$predicted - flu_bt$forecasts$observed),
                     na.rm = TRUE)
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
