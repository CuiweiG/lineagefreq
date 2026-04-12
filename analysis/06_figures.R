#!/usr/bin/env Rscript
# ============================================================
# 06_figures.R — Generate all publication figures
# ============================================================
# Loads results from 01-05 and produces all figures.
# Each figure uses cairo_pdf for correct rendering on Windows.
# All dimensions in millimeters, converted by save_figure().
# ============================================================

cat("[06_figures] Starting...\n")

if (!exists("theme_nature")) source("analysis/00_setup.R")

# ---- Load all results ----
bench <- readRDS("analysis/results/benchmark_ba2.rds")
cal   <- readRDS("analysis/results/calibration_results.rds")
surv  <- readRDS("analysis/results/surveillance_results.rds")
fitd  <- readRDS("analysis/results/fitness_results.rds")
flu   <- readRDS("analysis/results/influenza_results.rds")

detail <- bench$detail

# ============================================================
# Figure 1: Benchmark heatmap
# ============================================================
cat("[06_figures] Figure 1: Benchmark heatmap...\n")

fig1_data <- detail |>
  select(dataset, engine, horizon, mean_ae) |>
  mutate(label = sprintf("%.3f", mean_ae))

p1 <- ggplot(fig1_data, aes(x = factor(horizon), y = engine,
                             fill = mean_ae)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3, color = "white",
            fontface = "bold") +
  facet_wrap(~ dataset) +
  scale_fill_viridis(option = "C", direction = -1,
                     name = "Mean AE",
                     labels = scales::label_number(accuracy = 0.001)) +
  labs(x = "Forecast horizon (days)", y = "Engine",
       title = "Forecast accuracy: mean absolute error",
       caption = paste("Bedford Lab reference (2024):",
         "~0.6 pp median AE at 7d, ~6 pp mean AE at 14d.",
         if (length(bench$missing_engines) > 0)
           paste("\nNote:", paste(bench$missing_engines, collapse = ", "),
                 "unavailable (require CmdStan).")
         else ""))
save_figure(p1, "fig01_benchmark", width_mm = 180, height_mm = 120)


# ============================================================
# Figure 2: MAE by horizon with CI
# ============================================================
cat("[06_figures] Figure 2: MAE by horizon...\n")

ci_data <- detail |> filter(!is.na(mae_lower))

# Bedford Lab reference
bedford <- tibble(horizon = c(14, 28), bedford = c(0.06, 0.09))

p2 <- ggplot(ci_data, aes(x = horizon, y = mean_ae,
                           color = dataset, shape = engine)) +
  geom_pointrange(aes(ymin = mae_lower, ymax = mae_upper),
                  position = position_dodge(width = 2),
                  size = 0.5, linewidth = 0.4) +
  geom_hline(data = bedford, aes(yintercept = bedford),
             linetype = "dashed", color = "#999999", linewidth = 0.25) +
  annotate("text", x = 28.5, y = 0.06, label = "Bedford (14d)",
           size = 2.5, color = "#999999", hjust = 0) +
  annotate("text", x = 28.5, y = 0.09, label = "Bedford (28d)",
           size = 2.5, color = "#999999", hjust = 0) +
  scale_color_manual(values = c("BA.2" = "#0072B2",
                                "JN.1" = "#D55E00")) +
  scale_x_continuous(breaks = c(7, 14, 21, 28)) +
  labs(x = "Forecast horizon (days)",
       y = "Mean absolute error",
       color = "Dataset", shape = "Engine",
       title = "Forecast accuracy with 95% CI")
save_figure(p2, "fig02_mae_horizon", width_mm = 180, height_mm = 80)


# ============================================================
# Figure 3: PIT histograms
# ============================================================
cat("[06_figures] Figure 3: PIT histograms...\n")

pit_ba2 <- tibble(pit = cal$cal_ba2$pit_values, dataset = "BA.2")
pit_jn1 <- tibble(pit = cal$cal_jn1$pit_values, dataset = "JN.1")
pit_all <- bind_rows(pit_ba2, pit_jn1)

# Annotation data
ann <- tibble(
  dataset = c("BA.2", "JN.1"),
  label = c(
    sprintf("KS D = %.3f\np %s",
      cal$cal_ba2$ks_test$statistic,
      format.pval(cal$cal_ba2$ks_test$p_value, digits = 2)),
    sprintf("KS D = %.3f\np %s",
      cal$cal_jn1$ks_test$statistic,
      format.pval(cal$cal_jn1$ks_test$p_value, digits = 2))
  ),
  x = c(0.8, 0.8),
  y = c(Inf, Inf)
)

# Uniform reference line: if PIT values were U(0,1), each of 10
# bins would contain n/10 values. Compute per-dataset since
# facets use free_y and n differs.
n_ba2 <- length(cal$cal_ba2$pit_values)
n_jn1 <- length(cal$cal_jn1$pit_values)
ref_lines <- tibble(
  dataset  = c("BA.2", "JN.1"),
  yref     = c(n_ba2 / 10, n_jn1 / 10)
)

p3 <- ggplot(pit_all, aes(x = pit)) +
  geom_histogram(bins = 10, fill = "#6699CC", color = "white",
                 alpha = 0.85) +
  geom_hline(data = ref_lines, aes(yintercept = yref),
             linetype = "dashed", color = "#333333",
             linewidth = 0.25) +
  facet_wrap(~ dataset, scales = "free_y") +
  geom_text(data = ann, aes(x = x, y = y, label = label),
            vjust = 1.5, hjust = 0.5, size = 2.8,
            inherit.aes = FALSE) +
  labs(x = "PIT value", y = "Count",
       title = "PIT histograms: systematic underdispersion",
       subtitle = "U-shape indicates prediction intervals are too narrow")
save_figure(p3, "fig03_pit", width_mm = 180, height_mm = 80)


# ============================================================
# Figure 4: Reliability diagram
# ============================================================
cat("[06_figures] Figure 4: Reliability diagram...\n")

rel_all <- bind_rows(
  cal$cal_ba2$reliability |> mutate(dataset = "BA.2"),
  cal$cal_jn1$reliability |> mutate(dataset = "JN.1")
)

p4 <- ggplot(rel_all, aes(x = nominal, y = observed,
                           color = dataset)) +
  # Tolerance band
  geom_ribbon(data = tibble(nominal = seq(0, 1, 0.01)),
              aes(x = nominal, ymin = nominal - 0.10,
                  ymax = nominal + 0.10),
              fill = "#EEEEEE", color = NA, inherit.aes = FALSE) +
  # Perfect calibration
  geom_abline(slope = 1, intercept = 0, linetype = "solid",
              color = "#BBBBBB", linewidth = 0.25) +
  # Data
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  scale_color_manual(values = c("BA.2" = "#0072B2",
                                "JN.1" = "#D55E00")) +
  scale_x_continuous(limits = c(0, 1),
                     labels = scales::percent_format()) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format()) +
  coord_fixed() +
  labs(x = "Nominal coverage", y = "Observed coverage",
       color = "Dataset",
       title = "Reliability diagram",
       subtitle = "Gray band: +/- 10% tolerance")
save_figure(p4, "fig04_reliability", width_mm = 90, height_mm = 90)


# ============================================================
# Figure 5: Conformal vs parametric
# ============================================================
cat("[06_figures] Figure 5: Conformal vs parametric...\n")

# Select top 4 lineages by final frequency for clarity
data(cdc_sarscov2_jn1)
last_date <- max(cdc_sarscov2_jn1$date)
top_lins <- cdc_sarscov2_jn1 |>
  filter(date == last_date) |>
  arrange(desc(proportion)) |>
  head(4) |>
  pull(lineage)

param_sub <- cal$fc_parametric |>
  filter(.type == "forecast", .lineage %in% top_lins) |>
  mutate(method = "Parametric")
conf_sub <- cal$fc_conformal |>
  filter(.type == "forecast", .lineage %in% top_lins) |>
  mutate(method = "Conformal")
fc_cmp <- bind_rows(param_sub, conf_sub)

p5 <- ggplot(fc_cmp, aes(x = .date, y = .median,
                          color = .lineage)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper,
                  fill = .lineage),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ method) +
  scale_color_manual(values = pal_lineages[seq_along(top_lins)]) +
  scale_fill_manual(values = pal_lineages[seq_along(top_lins)]) +
  labs(x = "Date", y = "Frequency",
       color = "Lineage", fill = "Lineage",
       title = "Prediction intervals: parametric vs conformal",
       subtitle = "Conformal intervals are wider but correctly calibrated")
save_figure(p5, "fig05_conformal", width_mm = 180, height_mm = 100)


# ============================================================
# Figure 6: Surveillance optimization
# ============================================================
cat("[06_figures] Figure 6: Surveillance optimization...\n")

# Panel A: EVOI
p6a <- ggplot(surv$ev$values,
  aes(x = n_additional, y = evoi)) +
  geom_line(color = "#0072B2", linewidth = 0.5) +
  geom_point(color = "#0072B2", size = 1) +
  labs(x = "Additional sequences", y = "Variance reduction",
       title = "A. EVOI: diminishing returns")

# Panel B: Allocation
p6b <- plot(surv$ad_thompson) +
  labs(title = "B. Thompson sampling allocation")

p6 <- gridExtra::grid.arrange(p6a, p6b, ncol = 2)
save_figure(p6, "fig06_surveillance", width_mm = 180, height_mm = 70)


# ============================================================
# Figure 7: Alert threshold / JN.1 trajectory
# ============================================================
cat("[06_figures] Figure 7: Alert threshold...\n")

jn1_freq <- cdc_sarscov2_jn1[cdc_sarscov2_jn1$lineage == "JN.1", ]

p7 <- ggplot(jn1_freq, aes(x = date, y = proportion)) +
  geom_line(color = "#0072B2", linewidth = 0.6) +
  geom_point(color = "#0072B2", size = 1.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "#999999", linewidth = 0.25) +
  geom_vline(xintercept = surv$cross_5pct, linetype = "dotted",
             color = "#D55E00", linewidth = 0.3) +
  annotate("text", x = surv$cross_5pct + 14, y = 0.55,
           label = paste0("5% threshold\n", surv$cross_5pct),
           size = 2.5, color = "#D55E00", hjust = 0) +
  annotate("text", x = min(jn1_freq$date) + 10, y = 0.7,
           label = if (!surv$sprt_triggered)
             "SPRT: no alert on biweekly data\n(insufficient observations)"
           else paste0("SPRT alert: ", surv$alerts_sprt$date[
             surv$alerts_sprt$lineage == "JN.1"]),
           size = 2.5, color = "#333333", hjust = 0) +
  labs(x = "Date", y = "JN.1 proportion",
       title = "JN.1 emergence: surveillance detection",
       subtitle = "Biweekly CDC data")
save_figure(p7, "fig07_alert", width_mm = 180, height_mm = 70)


# ============================================================
# Figure 8: Fitness decomposition
# ============================================================
cat("[06_figures] Figure 8: Fitness decomposition...\n")

fd <- fitd$fd
d <- fd$decomposition[!is.na(fd$decomposition$transmissibility_fraction), ]

if (nrow(d) > 0) {
  plot_data <- d |>
    select(lineage, beta, escape_contribution) |>
    pivot_longer(cols = c(beta, escape_contribution),
                 names_to = "component", values_to = "value") |>
    mutate(component = ifelse(component == "beta",
      "Intrinsic transmissibility", "Immune escape"))

  p8 <- ggplot(plot_data, aes(x = lineage, y = value,
                               fill = component)) +
    geom_col(position = "stack", alpha = 0.85, width = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "#333333") +
    scale_fill_manual(values = c("Intrinsic transmissibility" = "#2166AC",
                                 "Immune escape" = "#B2182B")) +
    labs(x = "Lineage", y = "Growth rate component",
         fill = "Component",
         title = "Fitness decomposition: BA.2 transition",
         subtitle = "Generation time: 3.2 days (Du et al. 2022)",
         caption = "Positive escape: variant faces less population immunity than reference")
  save_figure(p8, "fig08_fitness", width_mm = 90, height_mm = 100)
} else {
  cat("  No decomposition data for non-pivot lineages.\n")
}


# ============================================================
# Figure 9: Influenza demonstration
# ============================================================
cat("[06_figures] Figure 9: Influenza demo...\n")

p9 <- autoplot(flu$fc_flu) +
  labs(title = "Influenza H3N2 clade dynamics",
       subtitle = "lineagefreq is pathogen-agnostic: same pipeline, different pathogen")
save_figure(p9, "fig09_influenza", width_mm = 180, height_mm = 70)


# ============================================================
# Figure 10: Graphical abstract
# ============================================================
cat("[06_figures] Figure 10: Graphical abstract...\n")

boxes <- tibble(
  x = 1:6,
  label = c("Data", "Model", "Calibrate",
            "Forecast", "Optimize", "Alert"),
  fn = c("lfq_data()", "fit_model()", "calibrate()",
         "forecast()", "adaptive_design()", "alert_threshold()")
)

p10 <- ggplot(boxes) +
  # Arrows
  geom_segment(data = boxes |> filter(x < 6),
    aes(x = x + 0.38, xend = x + 0.62, y = 0, yend = 0),
    arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
    color = "#666666", linewidth = 0.3) +
  # Boxes
  geom_label(aes(x = x, y = 0, label = label),
    fill = "#0072B2", color = "white", fontface = "bold",
    size = 3.2, label.padding = unit(3, "mm"),
    label.r = unit(1.5, "mm")) +
  # Function names
  geom_text(aes(x = x, y = -0.055, label = fn),
    size = 2.2, color = "#555555", family = "mono") +
  scale_x_continuous(limits = c(0.3, 6.7)) +
  scale_y_continuous(limits = c(-0.1, 0.08)) +
  labs(title = "lineagefreq: end-to-end genomic surveillance") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = 11))
save_figure(p10, "fig10_abstract", width_mm = 180, height_mm = 50)

cat("[06_figures] Complete.\n\n")
