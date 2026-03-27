# ============================================================
# Prepare built-in datasets for lineagefreq
# Run: Rscript data-raw/prepare_datasets.R
# ============================================================

# --- Dataset 1: Simulated SARS-CoV-2 US 2022 ---
# Realistic dynamics: BA.1 declining, BA.2 rising then declining,
# BA.4/5 emerging mid-year, BQ.1 emerging late.
# Simulated (not real GISAID data) to avoid license restrictions.

set.seed(20220101)

n_weeks  <- 40
dates    <- as.Date("2022-01-08") + (seq_len(n_weeks) - 1) * 7
variants <- c("BA.1", "BA.2", "BA.4/5", "BQ.1", "Other")
n_v      <- length(variants)

totals <- as.integer(pmin(
  pmax(round(rnorm(n_weeks, mean = 12000, sd = 4000)), 2000),
  25000
))

growth_rates <- c(
  "BA.1"   = -0.15,
  "BA.2"   =  0.08,
  "BA.4/5" =  0.25,
  "BQ.1"   =  0.20,
  "Other"  =  0.00
)

initial_logit <- c(
  "BA.1"   =  2.0,
  "BA.2"   = -1.0,
  "BA.4/5" = -6.0,
  "BQ.1"   = -8.0,
  "Other"  =  0.5
)

count_matrix <- matrix(NA_integer_, nrow = n_weeks, ncol = n_v)
for (t in seq_len(n_weeks)) {
  logits <- initial_logit + growth_rates * (t - 1)
  probs  <- exp(logits) / sum(exp(logits))
  count_matrix[t, ] <- stats::rmultinom(1, size = totals[t], prob = probs)
}

sarscov2_us_2022 <- data.frame(
  date    = rep(dates, each = n_v),
  variant = rep(variants, n_weeks),
  count   = as.vector(t(count_matrix)),
  total   = rep(totals, each = n_v),
  stringsAsFactors = FALSE
)

# --- Dataset 2: Simulated influenza A/H3N2 clade data ---

set.seed(20230901)

n_weeks_flu <- 24
dates_flu   <- as.Date("2023-10-01") + (seq_len(n_weeks_flu) - 1) * 7
clades      <- c("2a1b.2a.2", "2a1b.2a.1", "2a1b.1", "Other")
n_c         <- length(clades)

totals_flu <- as.integer(pmax(
  round(rnorm(n_weeks_flu, mean = 600, sd = 200)),
  100
))

growth_flu <- c(0.12, -0.05, -0.10, 0.00)
init_flu   <- c(-0.5,  1.0,   0.5,  0.0)

count_flu <- matrix(NA_integer_, nrow = n_weeks_flu, ncol = n_c)
for (t in seq_len(n_weeks_flu)) {
  logits <- init_flu + growth_flu * (t - 1)
  probs  <- exp(logits) / sum(exp(logits))
  count_flu[t, ] <- stats::rmultinom(1, size = totals_flu[t], prob = probs)
}

influenza_h3n2 <- data.frame(
  date  = rep(dates_flu, each = n_c),
  clade = rep(clades, n_weeks_flu),
  count = as.vector(t(count_flu)),
  total = rep(totals_flu, each = n_c),
  stringsAsFactors = FALSE
)

# --- Save ---
save(sarscov2_us_2022, file = "data/sarscov2_us_2022.rda", compress = "xz")
save(influenza_h3n2,   file = "data/influenza_h3n2.rda",   compress = "xz")

cat("Datasets created:\n")
cat(" sarscov2_us_2022:", nrow(sarscov2_us_2022), "rows\n")
cat(" influenza_h3n2:",   nrow(influenza_h3n2),   "rows\n")
