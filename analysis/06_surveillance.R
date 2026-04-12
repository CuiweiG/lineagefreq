###############################################################################
# 06_surveillance.R — Adaptive surveillance simulation (Conclusion 5)
# lineagefreq validation analysis
###############################################################################

cat("[06_surveillance] Starting...\n")
source("analysis/00_setup.R")

set.seed(20240416)

###############################################################################
# Section A: Adaptive allocation simulation
###############################################################################

cat("  [A] Adaptive allocation simulation...\n")

cdc_regional <- readRDS("analysis/results/cdc_regional_prepared.rds")
cat(sprintf("    Loaded %d HHS regional datasets\n", length(cdc_regional)))

# ── Simulation parameters ────────────────────────────────────────────────────
TOTAL_BUDGET   <- 2000L  # sequences per biweek nationally
N_REGIONS      <- 10L
N_REPLICATES   <- 50L
STRATEGIES     <- c("uniform", "proportional", "adaptive")

# HHS region approximate populations (millions, 2022 Census estimates)
# Used for proportional allocation
region_pop <- c(
  Region_1  = 15.1,   # CT, ME, MA, NH, RI, VT
  Region_2  = 22.2,   # NJ, NY
  Region_3  = 31.3,   # DE, DC, MD, PA, VA, WV
  Region_4  = 66.8,   # AL, FL, GA, KY, MS, NC, SC, TN
  Region_5  = 52.3,   # IL, IN, MI, MN, OH, WI
  Region_6  = 41.8,   # AR, LA, NM, OK, TX
  Region_7  = 14.3,   # IA, KS, MO, NE
  Region_8  = 12.4,   # CO, MT, ND, SD, UT, WY
  Region_9  = 52.7,   # AZ, CA, HI, NV
  Region_10 = 13.4    # AK, ID, OR, WA
)

# Normalize to match available regions
available_regions <- intersect(names(cdc_regional), names(region_pop))
region_pop <- region_pop[available_regions]
pop_total  <- sum(region_pop)

cat(sprintf("    %d regions with data, total budget = %d seq/biweek\n",
            length(available_regions), TOTAL_BUDGET))

# ── Extract "true" regional frequencies from CDC data ─────────────────────────

true_frequencies <- list()
common_dates     <- NULL

for (reg in available_regions) {
  obj <- cdc_regional[[reg]]
  props <- obj$counts / rowSums(obj$counts)
  true_frequencies[[reg]] <- list(
    dates = obj$dates,
    props = props,
    counts = obj$counts
  )
  if (is.null(common_dates)) {
    common_dates <- obj$dates
  } else {
    common_dates <- intersect(common_dates, obj$dates)
  }
}

common_dates <- sort(as.Date(common_dates, origin = "1970-01-01"))
n_timepoints <- length(common_dates)
cat(sprintf("    Common time points: %d\n", n_timepoints))

# ── Allocation functions ──────────────────────────────────────────────────────

allocate_uniform <- function(budget, n_regions) {
  rep(budget %/% n_regions, n_regions)
}

allocate_proportional <- function(budget, populations) {
  shares <- populations / sum(populations)
  alloc  <- round(shares * budget)
  # Adjust rounding to match budget exactly
  diff <- budget - sum(alloc)
  if (diff != 0) {
    idx <- order(shares * budget - alloc, decreasing = (diff > 0))
    alloc[idx[1:abs(diff)]] <- alloc[idx[1:abs(diff)]] + sign(diff)
  }
  alloc
}

allocate_adaptive <- function(budget, n_regions, entropy_scores) {
  # Thompson sampling-inspired: allocate proportional to information gain
  # Higher entropy (more lineage diversity / faster change) → more sequences
  # Softmax with temperature to avoid extreme allocations
  temperature <- 0.5
  weights     <- exp(entropy_scores / temperature)
  weights     <- weights / sum(weights)

  alloc <- round(weights * budget)
  diff  <- budget - sum(alloc)
  if (diff != 0) {
    idx <- order(weights * budget - alloc, decreasing = (diff > 0))
    alloc[idx[1:abs(diff)]] <- alloc[idx[1:abs(diff)]] + sign(diff)
  }
  pmax(alloc, 10L)  # minimum 10 sequences per region
}

# ── Simulation engine ─────────────────────────────────────────────────────────

run_simulation <- function(rep_id, strategy, budget, regions, true_freq,
                           common_dt, region_pops) {
  set.seed(20240416 + rep_id * 100 + match(strategy, STRATEGIES))

  n_reg <- length(regions)
  n_t   <- length(common_dt)

  mae_over_time   <- numeric(n_t)
  detection_times <- list()

  # Track entropy for adaptive allocation
  prev_entropy <- rep(1, n_reg)

  for (t_idx in seq_len(n_t)) {
    dt <- common_dt[t_idx]

    # Determine allocation
    alloc <- switch(strategy,
      uniform       = allocate_uniform(budget, n_reg),
      proportional  = allocate_proportional(budget, region_pops),
      adaptive      = allocate_adaptive(budget, n_reg, prev_entropy)
    )

    # Sample from true frequencies
    regional_errors <- numeric(n_reg)
    regional_preds  <- list()

    for (r_idx in seq_along(regions)) {
      reg  <- regions[r_idx]
      tf   <- true_freq[[reg]]

      # Find matching time point
      t_match <- which.min(abs(tf$dates - dt))
      true_props <- tf$props[t_match, ]
      true_props[is.na(true_props)] <- 0
      true_props <- true_props / max(sum(true_props), 1e-10)

      n_seq <- alloc[r_idx]
      n_lin <- length(true_props)

      # Multinomial sampling
      sampled_counts <- rmultinom(1, size = n_seq, prob = true_props)[, 1]
      sampled_props  <- sampled_counts / max(sum(sampled_counts), 1)

      regional_errors[r_idx] <- mean(abs(sampled_props - true_props))
      regional_preds[[reg]]  <- sampled_props

      # Update entropy for adaptive (Shannon entropy of sampled proportions)
      p <- sampled_props[sampled_props > 0]
      prev_entropy[r_idx] <- -sum(p * log(p))

      # Check detection: has BA.2 crossed 5% in this region?
      ba2_col <- grep("BA\\.2|BA2", colnames(tf$props), ignore.case = TRUE)
      if (length(ba2_col) > 0 && true_props[ba2_col[1]] > 0.05) {
        if (sampled_props[ba2_col[1]] > 0.05) {
          if (is.null(detection_times[[reg]])) {
            detection_times[[reg]] <- dt
          }
        }
      }
    }

    # National MAE: population-weighted average
    weights <- region_pops / sum(region_pops)
    mae_over_time[t_idx] <- sum(weights * regional_errors)
  }

  # Detection delay: compare to first date BA.2 truly > 5%
  detection_delays <- sapply(regions, function(reg) {
    tf <- true_freq[[reg]]
    ba2_col <- grep("BA\\.2|BA2", colnames(tf$props), ignore.case = TRUE)
    if (length(ba2_col) == 0) return(NA_real_)

    true_cross <- tf$dates[which(tf$props[, ba2_col[1]] > 0.05)[1]]
    detected   <- detection_times[[reg]]

    if (is.na(true_cross) || is.null(detected)) return(NA_real_)
    as.numeric(detected - true_cross)
  })

  tibble(
    replicate       = rep_id,
    strategy        = strategy,
    mae_trajectory  = list(mae_over_time),
    mean_mae        = mean(mae_over_time, na.rm = TRUE),
    detection_delay = list(detection_delays),
    mean_delay      = mean(detection_delays, na.rm = TRUE)
  )
}

# ── Run all simulations in parallel ───────────────────────────────────────────

sim_configs <- expand.grid(
  rep_id   = seq_len(N_REPLICATES),
  strategy = STRATEGIES,
  stringsAsFactors = FALSE
)

cat(sprintf("    Running %d simulations (%d strategies × %d reps)...\n",
            nrow(sim_configs), length(STRATEGIES), N_REPLICATES))

sim_results <- future_pmap_dfr(
  list(
    rep_id   = sim_configs$rep_id,
    strategy = sim_configs$strategy
  ),
  run_simulation,
  budget     = TOTAL_BUDGET,
  regions    = available_regions,
  true_freq  = true_frequencies,
  common_dt  = common_dates,
  region_pops = region_pop,
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

# Summarize by strategy
sim_summary <- sim_results |>
  group_by(strategy) |>
  summarise(
    mae_mean       = mean(mean_mae, na.rm = TRUE),
    mae_sd         = sd(mean_mae, na.rm = TRUE),
    delay_mean     = mean(mean_delay, na.rm = TRUE),
    delay_sd       = sd(mean_delay, na.rm = TRUE),
    n_reps         = n(),
    .groups = "drop"
  )

cat("    Simulation summary:\n")
print(sim_summary)

saveRDS(list(raw = sim_results, summary = sim_summary,
             config = list(budget = TOTAL_BUDGET, n_regions = N_REGIONS,
                           n_replicates = N_REPLICATES, strategies = STRATEGIES)),
        "analysis/results/surveillance_simulation.rds")
cat("    Saved surveillance_simulation.rds\n")

###############################################################################
# Section B: EVOI computation
###############################################################################

cat("  [B] EVOI computation...\n")

# Expected Value of Information: marginal value of additional sequencing
# Compute EVOI curve at different sample sizes

evoi_budgets <- c(100L, 200L, 500L, 1000L, 2000L)

# Use a representative BA.2 dataset
ba2_data <- tryCatch(cdc_sarscov2_ba2, error = function(e) NULL)

if (!is.null(ba2_data)) {
  evoi_results <- list()

  for (budget in evoi_budgets) {
    cat(sprintf("    EVOI at %d sequences...\n", budget))

    tryCatch({
      svd <- surveillance_value(ba2_data, n_sequences = budget)
      evoi_results[[as.character(budget)]] <- tibble(
        n_sequences = budget,
        evoi        = svd$evoi,
        mae_expected = svd$mae
      )
    },
    error = function(e) {
      cat(sprintf("      WARNING: surveillance_value() failed at n=%d: %s\n",
                  budget, e$message))
      # Fallback: estimate from sample size simulation results
      evoi_results[[as.character(budget)]] <<- tibble(
        n_sequences  = budget,
        evoi         = NA_real_,
        mae_expected = NA_real_
      )
    })
  }

  evoi_df <- bind_rows(evoi_results)

  saveRDS(evoi_df, "analysis/results/evoi_results.rds")
  cat("    Saved evoi_results.rds\n")
} else {
  cat("    Skipped EVOI (no BA.2 data)\n")
}

###############################################################################
# Section C: Alert threshold retrospective
###############################################################################

cat("  [C] Alert threshold retrospective...\n")

# Run alert_threshold() on JN.1 data to evaluate SPRT-based alerts

jn1_data <- tryCatch(cdc_sarscov2_jn1, error = function(e) NULL)

if (!is.null(jn1_data)) {
  tryCatch({
    alert_result <- alert_threshold(jn1_data)

    cat("    Alert threshold results:\n")
    print(alert_result)

    # Record alert date vs actual 5% crossing date
    jn1_props <- jn1_data$counts / rowSums(jn1_data$counts)
    jn1_cols  <- grep("JN\\.1|JN1", colnames(jn1_data$counts),
                       ignore.case = TRUE)

    if (length(jn1_cols) > 0) {
      actual_5pct <- jn1_data$dates[which(jn1_props[, jn1_cols[1]] > 0.05)[1]]
      cat(sprintf("    Actual 5%% crossing: %s\n", actual_5pct))

      if (!is.null(alert_result$alert_date)) {
        delay <- as.numeric(alert_result$alert_date - actual_5pct)
        cat(sprintf("    Alert date: %s (delay: %+d days)\n",
                    alert_result$alert_date, delay))
      }
    }

    saveRDS(alert_result, "analysis/results/alert_threshold.rds")
    cat("    Saved alert_threshold.rds\n")
  },
  error = function(e) {
    cat(sprintf("    WARNING: alert_threshold() failed: %s\n", e$message))
    cat("    This may be due to biweekly data granularity.\n")
    cat("    SPRT requires sufficient temporal resolution for sequential testing.\n")
    cat("    Consider interpolating to weekly data if available.\n")
  })
} else {
  cat("    Skipped alert threshold (no JN.1 data)\n")
}

cat("[06_surveillance] Complete.\n")
