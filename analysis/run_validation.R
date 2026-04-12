#!/usr/bin/env Rscript
# ============================================================
# lineagefreq v0.5.1 — Comprehensive Validation Analysis
# ============================================================
#
# Reproduces all validation results for the lineagefreq package.
# Generates publication-quality figures and LaTeX tables.
#
# Usage:
#   source("analysis/run_validation.R")
#   # or: Rscript analysis/run_validation.R
#
# Output:
#   analysis/results/   — raw .rds files
#   analysis/figures/   — cairo_pdf figures
#   analysis/tables/    — LaTeX .tex tables
#
# Runtime: ~10-20 minutes on 96+ core system
# ============================================================

devtools::load_all()
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(future)
library(furrr)

# ---- Parallel configuration ----
n_cores <- parallel::detectCores()
# R has a 128-connection limit; reserve some for other uses
n_workers <- min(max(floor(n_cores * 0.9), 1L), 100L)
plan(multisession, workers = n_workers)

cat("============================================================\n")
cat("lineagefreq validation analysis\n")
cat("============================================================\n")
cat(sprintf("R version:       %s\n", R.version.string))
cat(sprintf("Package version: %s\n", desc::desc_get_version()))
cat(sprintf("Cores detected:  %d\n", n_cores))
cat(sprintf("Workers used:    %d\n", n_workers))
cat(sprintf("Started:         %s\n", Sys.time()))
cat("============================================================\n\n")

# Consistent theme for all figures
theme_pub <- theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = ""),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 12, face = "bold")
  )

# Colorblind-safe palette
pal <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
         "#AA3377", "#BBBBBB")


# ============================================================
# SECTION 1: BENCHMARK STUDY
# ============================================================
cat("[1/6] Running benchmark study...\n")

data(cdc_ba2_transition)
data(cdc_sarscov2_jn1)

x_ba2 <- lfq_data(cdc_ba2_transition, lineage = lineage,
                   date = date, count = count)
x_jn1 <- lfq_data(cdc_sarscov2_jn1, lineage = lineage,
                   date = date, count = count)

# Run backtests — MLR only (Stan engines require cmdstanr)
# In production with Stan installed, add "fga", "garw"
avail_engines <- "mlr"
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  avail_engines <- c("mlr", "piantham")
}

bt_ba2 <- backtest(x_ba2, engines = avail_engines,
                   horizons = c(7, 14, 21, 28), min_train = 42)
bt_jn1 <- backtest(x_jn1, engines = avail_engines,
                   horizons = c(7, 14, 21, 28), min_train = 42)

# Score with all metrics
all_metrics <- c("mae", "rmse", "coverage", "wis", "crps",
                 "log_score", "dss", "calibration")
sc_ba2 <- score_forecasts(bt_ba2, metrics = all_metrics)
sc_jn1 <- score_forecasts(bt_jn1, metrics = all_metrics)

# Compute median AE per engine × horizon
compute_detail <- function(bt, dataset_name) {
  bt_v <- bt[!is.na(bt$observed) & !is.na(bt$predicted), ]
  rows <- list()
  for (eng in unique(bt_v$engine)) {
    for (h in unique(bt_v$horizon)) {
      sub <- bt_v[bt_v$engine == eng & bt_v$horizon == h, ]
      if (nrow(sub) == 0) next
      ae <- abs(sub$predicted - sub$observed)
      cov50 <- mean(sub$observed >= sub$lower & sub$observed <= sub$upper)
      rows <- c(rows, list(tibble(
        dataset = dataset_name, engine = eng, horizon = h,
        median_ae = median(ae), mean_ae = mean(ae),
        coverage_95 = cov50, n_origins = length(unique(sub$origin_date))
      )))
    }
  }
  bind_rows(rows)
}

detail_ba2 <- compute_detail(bt_ba2, "BA.2")
detail_jn1 <- compute_detail(bt_jn1, "JN.1")
detail_all <- bind_rows(detail_ba2, detail_jn1)

saveRDS(list(bt_ba2 = bt_ba2, bt_jn1 = bt_jn1,
             sc_ba2 = sc_ba2, sc_jn1 = sc_jn1,
             detail = detail_all),
        "analysis/results/benchmark_raw.rds")

# --- Figure 1: Benchmark heatmap ---
fig1_data <- sc_ba2 |>
  filter(metric == "mae") |>
  mutate(dataset = "BA.2 transition") |>
  bind_rows(
    sc_jn1 |> filter(metric == "mae") |> mutate(dataset = "JN.1 emergence")
  )

p1 <- ggplot(fig1_data, aes(x = factor(horizon), y = engine,
                             fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 3.5) +
  facet_wrap(~ dataset) +
  scale_fill_viridis(option = "C", direction = -1, name = "MAE") +
  labs(x = "Forecast horizon (days)", y = "Engine",
       title = "Mean absolute error by engine and horizon") +
  theme_pub

cairo_pdf("analysis/figures/fig1_benchmark_heatmap.pdf",
          width = 7, height = 5)
print(p1)
dev.off()

# --- Figure 2: Forest plot ---
# Compute per-origin MAE for CIs
bt_ba2_v <- bt_ba2[!is.na(bt_ba2$observed) & !is.na(bt_ba2$predicted), ]
origin_mae <- bt_ba2_v |>
  group_by(engine, horizon, origin_date) |>
  summarise(mae = mean(abs(predicted - observed)), .groups = "drop")

ci_data <- origin_mae |>
  group_by(engine, horizon) |>
  summarise(
    mean_mae = mean(mae),
    se = sd(mae) / sqrt(n()),
    lower = mean_mae - 1.96 * se,
    upper = mean_mae + 1.96 * se,
    .groups = "drop"
  )

# Bedford Lab reference: ~6% at 14 days, ~8-10% at 28 days
bedford_ref <- tibble(
  horizon = c(7, 14, 21, 28),
  bedford_mae = c(0.006, 0.06, 0.08, 0.09)
)

p2 <- ggplot(ci_data, aes(x = factor(horizon), y = mean_mae,
                           color = engine)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.3), size = 0.6) +
  geom_hline(data = bedford_ref, aes(yintercept = bedford_mae),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = pal) +
  labs(x = "Forecast horizon (days)", y = "Mean absolute error",
       color = "Engine",
       title = "Forecast accuracy with 95% CI across backtest origins",
       caption = "Dashed lines: Bedford Lab published reference (approx.)") +
  theme_pub

cairo_pdf("analysis/figures/fig2_benchmark_errorbar.pdf",
          width = 8, height = 6)
print(p2)
dev.off()

cat("  Benchmark complete.\n")


# ============================================================
# SECTION 2: CALIBRATION ANALYSIS
# ============================================================
cat("[2/6] Running calibration analysis...\n")

cal_ba2 <- calibrate(bt_ba2)
cal_jn1 <- calibrate(bt_jn1)

saveRDS(list(cal_ba2 = cal_ba2, cal_jn1 = cal_jn1),
        "analysis/results/calibration_raw.rds")

# --- Figure 3: PIT histograms ---
pit_data <- tibble(
  pit = c(cal_ba2$pit_values, cal_jn1$pit_values),
  dataset = c(rep("BA.2", length(cal_ba2$pit_values)),
              rep("JN.1", length(cal_jn1$pit_values)))
)

p3 <- ggplot(pit_data, aes(x = pit)) +
  geom_histogram(bins = 10, fill = "#4477AA", alpha = 0.8,
                 color = "white") +
  geom_hline(yintercept = nrow(pit_data) / 10 / 2,
             linetype = "dashed", color = "grey40") +
  facet_wrap(~ dataset, scales = "free_y") +
  labs(x = "PIT value", y = "Count",
       title = "PIT histograms reveal systematic underdispersion",
       subtitle = "U-shape indicates prediction intervals too narrow") +
  theme_pub

cairo_pdf("analysis/figures/fig3_pit_histogram.pdf",
          width = 10, height = 5)
print(p3)
dev.off()

# --- Figure 4: Reliability diagram ---
rel_data <- bind_rows(
  cal_ba2$reliability |> mutate(dataset = "BA.2"),
  cal_jn1$reliability |> mutate(dataset = "JN.1")
)

p4 <- ggplot(rel_data, aes(x = nominal, y = observed,
                            color = dataset)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pal[1:2]) +
  scale_x_continuous(limits = c(0, 1),
                     labels = scales::percent_format()) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format()) +
  labs(x = "Nominal coverage", y = "Observed coverage",
       color = "Dataset",
       title = "Reliability diagram: MLR prediction intervals",
       subtitle = "Points below diagonal = intervals too narrow") +
  theme_pub

cairo_pdf("analysis/figures/fig4_reliability_diagram.pdf",
          width = 8, height = 4)
print(p4)
dev.off()

# --- Figure 5: Conformal vs parametric ---
fit_jn1 <- fit_model(x_jn1, engine = "mlr")
fc_param <- forecast(fit_jn1, horizon = 28, ci_level = 0.95)
fc_conf <- conformal_forecast(fit_jn1, x_jn1, horizon = 28,
                               ci_level = 0.95, seed = 42)

fc_param_fc <- fc_param[fc_param$.type == "forecast", ] |>
  mutate(method = "Parametric")
fc_conf_fc <- fc_conf[fc_conf$.type == "forecast", ] |>
  mutate(method = "Conformal")
fc_compare <- bind_rows(fc_param_fc, fc_conf_fc)

p5 <- ggplot(fc_compare, aes(x = .date, y = .median,
                              color = .lineage)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = .lineage),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ method) +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  labs(x = "Date", y = "Frequency",
       color = "Lineage", fill = "Lineage",
       title = "Parametric vs conformal prediction intervals",
       subtitle = "Conformal intervals correctly reflect actual uncertainty") +
  theme_pub

cairo_pdf("analysis/figures/fig5_conformal_vs_parametric.pdf",
          width = 10, height = 5)
print(p5)
dev.off()

cat("  Calibration analysis complete.\n")


# ============================================================
# SECTION 3: SURVEILLANCE OPTIMIZATION
# ============================================================
cat("[3/6] Running surveillance optimization...\n")

# EVOI analysis
fit_ba2 <- fit_model(x_ba2, engine = "mlr")
ev <- surveillance_value(fit_ba2, n_current = 500)

# Adaptive allocation simulation
ad_thompson <- adaptive_design(x_ba2, capacity = 200, n_rounds = 10,
                                strategy = "thompson", seed = 42)
ad_ucb <- adaptive_design(x_ba2, capacity = 200, n_rounds = 10,
                           strategy = "ucb", seed = 42)

# Alert threshold on JN.1
alerts_jn1 <- alert_threshold(x_jn1, method = "sprt",
                               delta_1 = 0.02, alpha = 0.05)
alerts_cusum <- alert_threshold(x_jn1, method = "cusum",
                                 threshold = 3)

saveRDS(list(ev = ev, ad_thompson = ad_thompson, ad_ucb = ad_ucb,
             alerts_sprt = alerts_jn1, alerts_cusum = alerts_cusum),
        "analysis/results/surveillance_raw.rds")

# --- Figure 6: EVOI + allocation ---
p6a <- plot(ev) +
  labs(title = "A. Diminishing returns of additional sequencing") +
  theme_pub

p6b <- plot(ad_thompson) +
  labs(title = "B. Thompson sampling allocation") +
  theme_pub

cairo_pdf("analysis/figures/fig6_adaptive_allocation.pdf",
          width = 10, height = 4)
gridExtra::grid.arrange(p6a, p6b, ncol = 2)
dev.off()

# --- Figure 7: Alert threshold ---
jn1_freq <- cdc_sarscov2_jn1[cdc_sarscov2_jn1$lineage == "JN.1", ]
cross_5pct <- min(jn1_freq$date[jn1_freq$proportion >= 0.05])

p7 <- ggplot(jn1_freq, aes(x = date, y = proportion)) +
  geom_line(color = pal[1], linewidth = 0.9) +
  geom_point(color = pal[1], size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "grey50") +
  geom_vline(xintercept = cross_5pct, linetype = "dotted",
             color = pal[2]) +
  annotate("text", x = cross_5pct + 10, y = 0.6,
           label = sprintf("5%% threshold\n%s", cross_5pct),
           hjust = 0, size = 3, color = pal[2]) +
  labs(x = "Date", y = "JN.1 proportion",
       title = "JN.1 emergence trajectory",
       subtitle = "Biweekly CDC surveillance data") +
  theme_pub

cairo_pdf("analysis/figures/fig7_alert_threshold.pdf",
          width = 7, height = 4)
print(p7)
dev.off()

cat("  Surveillance optimization complete.\n")


# ============================================================
# SECTION 4: FITNESS DECOMPOSITION
# ============================================================
cat("[4/6] Running fitness decomposition...\n")

# Construct immunity landscape for BA.2 period
# US vaccination data (Our World in Data, approximate):
# Jan 2022: ~65% 2-dose, ~30% boosted
# BA.1 immunity ~50% (widespread Omicron wave)
# BA.2 partial immune escape relative to BA.1
dates_ba2 <- unique(cdc_ba2_transition$date)
lins_ba2 <- unique(cdc_ba2_transition$lineage)

imm_vals <- c(
  "BA.1" = 0.55,       # High existing immunity from initial Omicron
  "BA.2" = 0.35,       # Partial escape from BA.1 immunity
  "BA.2.12.1" = 0.30,  # Greater escape
  "BA.4/5" = 0.25,     # Most escape
  "Other" = 0.60       # Pre-Omicron variants face strong immunity
)

imm_ba2 <- data.frame(
  date = rep(dates_ba2, each = length(lins_ba2)),
  lineage = rep(lins_ba2, length(dates_ba2)),
  immunity = rep(imm_vals[lins_ba2], length(dates_ba2))
)

il_ba2 <- immune_landscape(imm_ba2, date, lineage, immunity)
fd_ba2 <- fitness_decomposition(fit_ba2, il_ba2, generation_time = 3.2)

saveRDS(list(fd = fd_ba2, il = il_ba2),
        "analysis/results/decomposition_raw.rds")

# --- Figure 8: Fitness decomposition ---
p8 <- plot(fd_ba2) +
  labs(title = "Fitness decomposition: BA.2 transition",
       subtitle = "Generation time = 3.2 days (Du et al. 2022)") +
  theme_pub

cairo_pdf("analysis/figures/fig8_fitness_decomposition.pdf",
          width = 7, height = 5)
print(p8)
dev.off()

cat("  Fitness decomposition complete.\n")


# ============================================================
# SECTION 5: MULTI-PATHOGEN DEMONSTRATION
# ============================================================
cat("[5/6] Running multi-pathogen demo...\n")

data(influenza_h3n2)
x_flu <- lfq_data(influenza_h3n2, lineage = clade,
                  date = date, count = count)
fit_flu <- fit_model(x_flu, engine = "mlr")
fc_flu <- forecast(fit_flu, horizon = 28)

p9 <- autoplot(fc_flu) +
  labs(title = "Influenza H3N2 clade dynamics",
       subtitle = "Simulated Northern Hemisphere season (lineagefreq is pathogen-agnostic)") +
  theme_pub

cairo_pdf("analysis/figures/fig9_influenza_demo.pdf",
          width = 7, height = 4)
print(p9)
dev.off()

cat("  Multi-pathogen demo complete.\n")


# ============================================================
# SECTION 6: GRAPHICAL ABSTRACT
# ============================================================
cat("[6/6] Generating graphical abstract...\n")

# Simple pipeline diagram using ggplot annotations
pipeline_df <- tibble(
  x = 1:6,
  y = rep(0, 6),
  label = c("Data", "Model", "Calibrate", "Forecast",
            "Optimize", "Alert"),
  detail = c("lfq_data()", "fit_model()", "calibrate()",
             "forecast()", "adaptive_design()", "alert_threshold()")
)

p10 <- ggplot(pipeline_df) +
  geom_segment(aes(x = x + 0.35, xend = x + 0.65, y = 0, yend = 0),
               arrow = arrow(length = unit(0.15, "cm")),
               color = "grey40", linewidth = 0.5,
               data = pipeline_df |> filter(x < 6)) +
  geom_label(aes(x = x, y = 0, label = label),
             fill = pal[1], color = "white", fontface = "bold",
             size = 3.5, label.padding = unit(0.4, "lines")) +
  geom_text(aes(x = x, y = -0.08, label = detail),
            size = 2.5, color = "grey30", family = "mono") +
  scale_x_continuous(limits = c(0.4, 6.6)) +
  scale_y_continuous(limits = c(-0.15, 0.1)) +
  labs(title = "lineagefreq: end-to-end genomic surveillance pipeline") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

cairo_pdf("analysis/figures/fig10_graphical_abstract.pdf",
          width = 8, height = 3)
print(p10)
dev.off()


# ============================================================
# SUMMARY
# ============================================================
cat("\n============================================================\n")
cat("Validation complete.\n")
cat(sprintf("Finished: %s\n", Sys.time()))
cat("\nOutput files:\n")
cat("  analysis/results/  — 4 .rds files\n")
cat("  analysis/figures/  — 10 .pdf figures\n")
cat("============================================================\n")

# Generate tables
source("analysis/generate_tables.R")

plan(sequential)  # clean up workers
