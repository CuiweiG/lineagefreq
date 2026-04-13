###############################################################################
# 13_caps_theorem.R
# Empirical verification of CAPS coverage guarantee
#
# Run from package root:
#   source("analysis/13_caps_theorem.R")
#
# Requires: analysis/results/caps_comparison.rds (from 10_caps_validation.R)
###############################################################################

cat("[13_caps_theorem] Starting...\n")
source("analysis/00_setup.R")
library(lineagefreq)

# ─── Theorem statement (verified empirically) ────────────────────────────
#
# Let {(X_t, Y_t)} be a non-stationary time series with Y_t in Delta^{K-1}
# (the K-simplex). Let R_hat(h) be the variance ratio estimated from
# n_cal calibration samples. CAPS with radius
#
#   r(t, h) = r_base(h) * psi(K)
#
# satisfies the coverage bound:
#
#   P(Y_{t+1}(h) in C_CAPS(h)) >= 1 - alpha - eps_finite - eps_shift
#
# where:
#   eps_finite(n_cal) = (1 - alpha) / n_cal    [finite-sample correction]
#   eps_shift = gamma(h) * sup_t |alpha_t - alpha|  [ACI regret]
#
# Key insight: gamma(h) = gamma_0 * max(0, 1 - R_hat(h))
# When R_hat(h) ~ 1 (short horizon): gamma ~ 0, eps_shift ~ 0
# When R_hat(h) << 1 (long horizon): gamma > 0 but r_base already
# compensates via the sqrt(1/R) factor, so net coverage stays near 1-alpha.

# ─── Load CAPS results ──────────────────────────────────────────────────

caps_file <- "analysis/results/caps_comparison.rds"
if (!file.exists(caps_file)) {
  cat("  caps_comparison.rds not found. Run 10_caps_validation.R first.\n")
  cat("[13_caps_theorem] Skipped.\n")
  quit(save = "no")
}

caps <- readRDS(caps_file)

if (is.null(caps) || nrow(caps) == 0) {
  cat("  No CAPS comparison data available.\n")
  cat("[13_caps_theorem] Skipped.\n")
  quit(save = "no")
}

cat(sprintf("  Loaded %d dataset-horizon pairs\n", nrow(caps)))

# ─── Compute theoretical coverage bounds ─────────────────────────────────

alpha <- 0.05
gamma_0 <- 0.05

caps$theoretical_bound <- NA_real_
caps$eps_finite <- NA_real_
caps$eps_shift <- NA_real_
caps$gamma_h <- NA_real_
caps$bound_satisfied <- NA

cat("\n-- CAPS coverage guarantee verification --\n\n")

for (i in seq_len(nrow(caps))) {
  n_cal <- caps$n_test[i]  # approximate calibration set size
  R_h   <- caps$R_hat[i]

  if (is.na(R_h) || is.na(n_cal) || n_cal < 2) next

  gamma_h <- gamma_0 * max(0, 1 - R_h)

  # Finite-sample correction
  eps_finite <- (1 - alpha) / n_cal

  # ACI regret bound (conservative: sup |alpha_t - alpha| <= 0.1)
  eps_shift <- gamma_h * 0.1

  theoretical_lower <- 1 - alpha - eps_finite - eps_shift

  caps$eps_finite[i]        <- eps_finite
  caps$eps_shift[i]         <- eps_shift
  caps$gamma_h[i]           <- gamma_h
  caps$theoretical_bound[i] <- theoretical_lower
  caps$bound_satisfied[i]   <- caps$caps_cov[i] >= theoretical_lower

  cat(sprintf("  h=%2dd: R=%.2f gamma=%.3f n=%3d | ",
              caps$horizon[i], R_h, gamma_h, n_cal))
  cat(sprintf("Observed=%.1f%% >= Bound=%.1f%% ? %s\n",
              caps$caps_cov[i] * 100, theoretical_lower * 100,
              ifelse(caps$bound_satisfied[i], "YES", "NO")))
}

# ─── Summary ─────────────────────────────────────────────────────────────

valid_rows <- !is.na(caps$bound_satisfied)
n_valid    <- sum(valid_rows)
n_satisfied <- sum(caps$bound_satisfied, na.rm = TRUE)

cat(sprintf("\n-- Summary --\n"))
cat(sprintf("Valid comparisons: %d\n", n_valid))
cat(sprintf("Theoretical bound satisfied: %d/%d (%.0f%%)\n",
            n_satisfied, n_valid, n_satisfied / n_valid * 100))
cat(sprintf("Mean observed CAPS coverage: %.1f%%\n",
            mean(caps$caps_cov[valid_rows], na.rm = TRUE) * 100))
cat(sprintf("Mean theoretical lower bound: %.1f%%\n",
            mean(caps$theoretical_bound[valid_rows], na.rm = TRUE) * 100))
cat(sprintf("Excess coverage (observed - bound): %.1f pp\n",
            mean(caps$caps_cov[valid_rows] - caps$theoretical_bound[valid_rows],
                 na.rm = TRUE) * 100))

# ─── Horizon-dependent analysis ──────────────────────────────────────────

cat("\n-- Coverage bound by horizon --\n")
for (h in sort(unique(caps$horizon))) {
  h_rows <- caps$horizon == h & valid_rows
  if (sum(h_rows) == 0) next

  obs_mean  <- mean(caps$caps_cov[h_rows], na.rm = TRUE)
  bound_mean <- mean(caps$theoretical_bound[h_rows], na.rm = TRUE)
  sat_frac  <- mean(caps$bound_satisfied[h_rows], na.rm = TRUE)

  cat(sprintf("  h=%2dd: observed=%.1f%%, bound=%.1f%%, satisfied=%d/%d\n",
              h, obs_mean * 100, bound_mean * 100,
              sum(caps$bound_satisfied[h_rows], na.rm = TRUE), sum(h_rows)))
}

# ─── Save ────────────────────────────────────────────────────────────────

saveRDS(caps, "analysis/results/caps_theorem_verification.rds")
cat("\nSaved analysis/results/caps_theorem_verification.rds\n")
cat("\n[13_caps_theorem] Complete.\n")
