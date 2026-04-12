###############################################################################
# 08_figures.R — Publication figures (5 focused figures + supplementary)
# lineagefreq validation analysis
# Restructured per expert review: each figure tells ONE story, 2-4 panels max
###############################################################################

cat("[08_figures] Starting...\n")
source("analysis/00_setup.R")

# ─── Load all results ────────────────────────────────────────────────────────

benchmark   <- readRDS("analysis/results/benchmark_multicountry.rds")
calibration <- readRDS("analysis/results/calibration_comparison.rds")
decision    <- readRDS("analysis/results/decision_impact.rds")
sample_size <- tryCatch(readRDS("analysis/results/sample_size_analysis.rds"), error = function(e) NULL)
window_anal <- tryCatch(readRDS("analysis/results/window_analysis.rds"), error = function(e) NULL)
fitness     <- tryCatch(readRDS("analysis/results/fitness_sensitivity.rds"), error = function(e) NULL)
surveillance <- readRDS("analysis/results/surveillance_simulation.rds")
evoi        <- tryCatch(readRDS("analysis/results/evoi_results.rds"), error = function(e) NULL)
influenza   <- readRDS("analysis/results/influenza_results.rds")

# Built-in data: data.frames with columns date, lineage, count, proportion
ba2_data <- tryCatch(cdc_ba2_transition, error = function(e) NULL)
jn1_data <- tryCatch(cdc_sarscov2_jn1, error = function(e) NULL)

# ─── Okabe-Ito colorblind-safe palette ───────────────────────────────────────
# Never use red+green together
oi <- c(
  orange    = "#E69F00",
  sky_blue  = "#56B4E9",
  green     = "#009E73",
  yellow    = "#F0E442",
  dark_blue = "#0072B2",
  vermillion = "#D55E00",
  pink      = "#CC79A7",
  gray      = "#999999"
)

# Extended Okabe-Ito palette — enough for up to 10 categories
pal_extended <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9",
                  "#CC79A7", "#F0E442", "#999999", "#000000", "#332288")

# Dynamic palette: returns n colours from the extended set
pal_n <- function(n) head(pal_extended, max(n, 1))

# 2-method comparison: Parametric vs Conformal (NO Recalibrated)
pal_2methods <- c(Parametric = unname(oi["dark_blue"]),
                  Conformal  = unname(oi["vermillion"]))

# Fitness components
pal_fitness <- c(Transmissibility = unname(oi["dark_blue"]),
                 Immune_escape    = unname(oi["orange"]))

# ─── Consistent panel label helper ──────────────────────────────────────────
# Bold 10pt, top-left corner, fixed position
panel_label <- function(label) {
  annotate("text", x = -Inf, y = Inf, label = label,
           fontface = "bold", family = "Arial", size = 10 / .pt,
           hjust = -0.3, vjust = 1.5)
}

###############################################################################
# FIGURE 1: Method and Core Capability (180 x 120 mm)
###############################################################################

cat("  [Figure 1] Method and core capability...\n")

# ── Panel a: BA.2 frequency dynamics ─────────────────────────────────────────

fig1a <- tryCatch({
  if (!is.null(ba2_data)) {
    # Compute frequencies, keep top 5 lineages by max proportion
    freq_df <- ba2_data |>
      group_by(date, lineage) |>
      summarise(count = sum(count), .groups = "drop") |>
      group_by(date) |>
      mutate(frequency = count / sum(count)) |>
      ungroup()

    top5 <- freq_df |>
      group_by(lineage) |>
      summarise(peak = max(frequency, na.rm = TRUE), .groups = "drop") |>
      slice_max(peak, n = 5) |>
      pull(lineage)

    plot_df <- freq_df |> filter(lineage %in% top5)

    ggplot(plot_df, aes(x = date, y = frequency, colour = lineage)) +
      geom_point(size = 0.8, alpha = 0.5) +
      geom_line(linewidth = LW_DATA) +
      scale_colour_manual(values = setNames(pal_n(length(top5)), top5)) +
      scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 month") +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      labs(x = "Date", y = "Frequency", colour = NULL) +
      panel_label("a") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            legend.key.width = unit(5, "mm"))
  } else {
    ggplot() + annotate("text", x = .5, y = .5, label = "Data unavailable") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel b: Pipeline schematic ──────────────────────────────────────────────

fig1b <- ggplot() +
  xlim(0, 12) + ylim(0, 4) +
  # Box 1: Data
  annotate("rect", xmin = 0.2, xmax = 2.3, ymin = 1.5, ymax = 2.5,
           fill = "#E8F0FE", colour = "black", linewidth = 0.3) +
  annotate("text", x = 1.25, y = 2.1, label = "Data",
           family = "Arial", size = 8 / .pt, fontface = "bold") +
  annotate("text", x = 1.25, y = 1.7, label = "lfq_data()",
           family = "Arial", size = 5.5 / .pt, colour = "grey40") +
  # Arrow 1-2
  annotate("segment", x = 2.3, xend = 3.0, y = 2.0, yend = 2.0,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.3) +
  # Box 2: Model
  annotate("rect", xmin = 3.0, xmax = 5.2, ymin = 1.5, ymax = 2.5,
           fill = "#E8F0FE", colour = "black", linewidth = 0.3) +
  annotate("text", x = 4.1, y = 2.1, label = "Model",
           family = "Arial", size = 8 / .pt, fontface = "bold") +
  annotate("text", x = 4.1, y = 1.7, label = "fit_model()",
           family = "Arial", size = 5.5 / .pt, colour = "grey40") +
  # Arrow 2-3
  annotate("segment", x = 5.2, xend = 5.9, y = 2.0, yend = 2.0,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.3) +
  # Box 3: Forecast
  annotate("rect", xmin = 5.9, xmax = 8.1, ymin = 1.5, ymax = 2.5,
           fill = "#E8F0FE", colour = "black", linewidth = 0.3) +
  annotate("text", x = 7.0, y = 2.1, label = "Forecast",
           family = "Arial", size = 8 / .pt, fontface = "bold") +
  annotate("text", x = 7.0, y = 1.7, label = "forecast()",
           family = "Arial", size = 5.5 / .pt, colour = "grey40") +
  # Arrow 3-4
  annotate("segment", x = 8.1, xend = 8.8, y = 2.0, yend = 2.0,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.3) +
  # Box 4: Calibrate (dashed border = optional step)
  annotate("rect", xmin = 8.8, xmax = 11.2, ymin = 1.5, ymax = 2.5,
           fill = "#F5F0FF", colour = "black", linewidth = 0.3,
           linetype = "dashed") +
  annotate("text", x = 10.0, y = 2.1, label = "Calibrate",
           family = "Arial", size = 8 / .pt, fontface = "bold") +
  annotate("text", x = 10.0, y = 1.7, label = "calibrate()",
           family = "Arial", size = 5.5 / .pt, colour = "grey40") +
  panel_label("b") +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

figure1 <- fig1a + fig1b + plot_layout(widths = c(3, 2))
save_figure(figure1, "figure1", width_mm = 180, height_mm = 120)

###############################################################################
# FIGURE 2: Multi-country Validation (180 x 100 mm)
###############################################################################

cat("  [Figure 2] Multi-country validation...\n")

# ── Panel a: 5-country frequency dynamics ────────────────────────────────────

fig2a <- tryCatch({
  ecdc_data <- readRDS("analysis/results/ecdc_prepared.rds")

  ecdc_plot_list <- list()
  for (nm in names(ecdc_data)) {
    if (!grepl("BA2", nm)) next
    obj <- ecdc_data[[nm]]
    country <- sub("_BA2$", "", nm)

    # lfq_data tibble: .date, .lineage, .freq
    cdf <- obj |>
      select(date = .date, lineage = .lineage, frequency = .freq) |>
      mutate(country = country)

    # Keep top 4 lineages per country
    top4 <- cdf |>
      group_by(lineage) |>
      summarise(peak = max(frequency, na.rm = TRUE), .groups = "drop") |>
      slice_max(peak, n = 4) |>
      pull(lineage)
    ecdc_plot_list[[nm]] <- cdf |> filter(lineage %in% top4)
  }

  if (length(ecdc_plot_list) > 0) {
    ecdc_df <- bind_rows(ecdc_plot_list)

    n_lin <- n_distinct(ecdc_df$lineage)
    ggplot(ecdc_df, aes(x = date, y = frequency, colour = lineage)) +
      geom_line(linewidth = LW_DATA) +
      facet_wrap(~country, nrow = 1) +
      scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "3 months") +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      scale_colour_manual(values = pal_n(n_lin)) +
      labs(x = "Date", y = "Frequency", colour = NULL) +
      panel_label("a") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            strip.text = element_text(size = 7),
            legend.position = "bottom",
            legend.key.width = unit(4, "mm"))
  } else {
    ggplot() + annotate("text", x = .5, y = .5, label = "ECDC data unavailable") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel b: MAE benchmark forest plot ───────────────────────────────────────

fig2b <- tryCatch({
  mae_df <- benchmark$metrics |>
    filter(engine == "mlr", horizon == 14) |>
    left_join(benchmark$bootstrap_ci |> filter(engine == "mlr", horizon == 14),
              by = c("dataset", "engine", "horizon")) |>
    mutate(mae_pct = mae * 100,
           lo_pct  = mae_lower * 100,
           hi_pct  = mae_upper * 100,
           model   = "MLR")

  # Compute naive baseline MAE: last observed frequency carried forward
  # For each backtest origin, naive_predicted = observed at origin_date
  naive_mae <- map_dfr(names(benchmark$backtest_results), function(nm) {
    bt <- benchmark$backtest_results[[nm]]
    if (bt$engine != "mlr") return(NULL)
    bt_tbl <- bt$result |> filter(horizon == 14)
    if (nrow(bt_tbl) == 0) return(NULL)

    # For each origin_date × lineage, find the observed frequency at origin_date
    # (the most recent observation before the forecast)
    # Naive prediction = predicted freq at horizon 0 = the last known observed
    # We approximate: naive error ≈ |observed_at_target - observed_at_origin|
    # Since we don't have obs at origin directly, use predicted as proxy for
    # a well-fitted model; for the naive model, the error is typically larger.
    # Better: compute from all origins — the last-value-carried-forward MAE
    naive_errors <- bt_tbl |>
      group_by(origin_date, lineage) |>
      slice_head(n = 1) |>
      ungroup()

    # Get observed at origin_date from a different horizon row or from data
    # Simplest valid approach: for each (origin, lineage) pair, the naive
    # forecast at target_date is the observed at origin_date. We can get this
    # from horizon==7 rows where target_date is closer to origin_date.
    bt_all <- bt$result
    origin_obs <- bt_all |>
      filter(horizon == min(horizon)) |>
      # The observed at the closest target is a proxy for frequency at origin
      select(origin_date, lineage, naive_pred = observed)

    naive_joined <- bt_tbl |>
      left_join(origin_obs, by = c("origin_date", "lineage")) |>
      filter(!is.na(naive_pred), !is.na(observed))

    if (nrow(naive_joined) == 0) return(NULL)

    tibble(
      dataset = bt$dataset,
      mae_pct = mean(abs(naive_joined$naive_pred - naive_joined$observed)) * 100
    )
  })

  # Combine MLR and naive into one plot data frame
  plot_df <- mae_df |>
    select(dataset, mae_pct, lo_pct, hi_pct, model)

  if (nrow(naive_mae) > 0) {
    naive_plot <- naive_mae |>
      mutate(model = "Naive", lo_pct = NA_real_, hi_pct = NA_real_)
    plot_df <- bind_rows(plot_df, naive_plot)
  }

  ggplot(plot_df, aes(x = mae_pct, y = reorder(dataset, -mae_pct),
                      colour = model, shape = model)) +
    geom_point(size = 2, position = position_dodge(width = 0.4)) +
    geom_errorbarh(aes(xmin = lo_pct, xmax = hi_pct), height = 0.3,
                   linewidth = LW_DATA,
                   position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = c(MLR = oi["dark_blue"],
                                   Naive = oi["gray"])) +
    scale_shape_manual(values = c(MLR = 16, Naive = 17)) +
    labs(x = "MAE at 14-day horizon (%)", y = NULL, colour = NULL, shape = NULL) +
    panel_label("b") +
    theme_nature() +
    theme(legend.position = "bottom",
          legend.key.size = unit(3, "mm"))
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

figure2 <- fig2a + fig2b + plot_layout(widths = c(3, 2))
save_figure(figure2, "figure2", width_mm = 180, height_mm = 100)

###############################################################################
# FIGURE 3: Calibration — The Core Finding (180 x 140 mm)
###############################################################################

cat("  [Figure 3] Calibration...\n")

# ── Panel a: PIT histogram (US BA.2 only) ────────────────────────────────────

fig3a <- tryCatch({
  # Find US BA.2 mlr PIT result
  ba2_pit_key <- names(calibration$pit_results) |>
    grep("BA2.*mlr|ba2.*mlr|BA2_builtin.*mlr", x = _, value = TRUE)
  if (length(ba2_pit_key) == 0) ba2_pit_key <- names(calibration$pit_results)[1]

  pr <- calibration$pit_results[[ba2_pit_key[1]]]

  ggplot(tibble(pit = pr$pit_values), aes(x = pit)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = oi["dark_blue"], colour = "white", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = LW_REF,
               colour = oi["gray"]) +
    annotate("text", x = 0.95, y = Inf,
             label = sprintf("KS D = %.3f\np = %.2g", pr$ks_D, pr$ks_p),
             hjust = 1, vjust = 1.3, size = 6 / .pt, family = "Arial") +
    labs(x = "PIT value", y = "Density") +
    panel_label("a") +
    theme_nature()
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel b: Reliability diagram ─────────────────────────────────────────────

fig3b <- tryCatch({
  rel_df <- calibration$reliability |>
    filter(engine == "mlr", !is.na(observed_coverage))

  ggplot(rel_df, aes(x = nominal_coverage, y = observed_coverage,
                     colour = dataset, group = dataset)) +
    # Tolerance band
    geom_ribbon(data = tibble(x = seq(0, 1, 0.01)),
                aes(x = x, ymin = pmax(0, x - 0.10), ymax = pmin(1, x + 0.10)),
                inherit.aes = FALSE, fill = "grey85", alpha = 0.4) +
    # Diagonal
    geom_abline(slope = 1, intercept = 0, linewidth = LW_REF, colour = "black") +
    geom_line(linewidth = LW_DATA) +
    geom_point(size = 1) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_colour_manual(values = pal_n(n_distinct(rel_df$dataset))) +
    labs(x = "Nominal coverage", y = "Observed coverage", colour = NULL) +
    panel_label("b") +
    theme_nature() +
    theme(legend.position = "bottom",
          legend.key.width = unit(4, "mm"))
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel c: Parametric vs Conformal intervals ──────────────────────────────
# TWO sub-rows only: Parametric (top) and Conformal (bottom)
# NO Recalibrated — isotonic regression produces numerical overflow

fig3c <- tryCatch({
  # Find first available calibration comparison (prefer BA2)
  cc_keys <- names(calibration$calibration_comparison)
  ba2_cc <- grep("BA2|ba2", cc_keys, value = TRUE)
  cc_key <- if (length(ba2_cc) > 0) ba2_cc[1] else cc_keys[1]
  cal_comp <- calibration$calibration_comparison[[cc_key]]

  if (!is.null(cal_comp$parametric) && !is.null(cal_comp$conformal)) {
    para_df <- cal_comp$parametric
    conf_df <- cal_comp$conformal

    # Pick dominant lineage, shortest horizon
    top_lin <- para_df |>
      group_by(lineage) |>
      summarise(peak = max(observed, na.rm = TRUE), .groups = "drop") |>
      slice_max(peak, n = 1) |>
      pull(lineage)

    min_h <- min(para_df$horizon, na.rm = TRUE)

    para_lin <- para_df |> filter(lineage == top_lin[1], horizon == min_h)
    conf_lin <- conf_df |> filter(lineage == top_lin[1], horizon == min_h)

    # Parametric: observations outside lower/upper
    para_outside <- para_lin |> filter(observed < lower | observed > upper)

    # Conformal: observations outside conf_lower/conf_upper
    conf_outside <- if ("conf_lower" %in% names(conf_lin)) {
      conf_lin |> filter(observed < conf_lower | observed > conf_upper)
    } else tibble()

    p_para <- ggplot(para_lin, aes(x = origin_date)) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  fill = pal_2methods["Parametric"], alpha = 0.2) +
      geom_line(aes(y = predicted), colour = pal_2methods["Parametric"],
                linewidth = LW_DATA) +
      geom_point(aes(y = observed), size = 0.5, colour = "black") +
      {if (nrow(para_outside) > 0)
        geom_point(data = para_outside, aes(y = observed),
                   size = 0.8, colour = oi["vermillion"], shape = 16)} +
      scale_x_date(date_labels = "%Y-%m-%d") +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = NULL, y = "Frequency", title = "Parametric") +
      theme_nature() +
      theme(axis.text.x = element_blank(),
            plot.title = element_text(size = 8))

    lo_col <- if ("conf_lower" %in% names(conf_lin)) "conf_lower" else "lower"
    hi_col <- if ("conf_upper" %in% names(conf_lin)) "conf_upper" else "upper"

    p_conf <- ggplot(conf_lin, aes(x = origin_date)) +
      geom_ribbon(aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]),
                  fill = pal_2methods["Conformal"], alpha = 0.2) +
      geom_line(aes(y = predicted), colour = pal_2methods["Conformal"],
                linewidth = LW_DATA) +
      geom_point(aes(y = observed), size = 0.5, colour = "black") +
      {if (nrow(conf_outside) > 0)
        geom_point(data = conf_outside, aes(y = observed),
                   size = 0.8, colour = oi["vermillion"], shape = 16)} +
      scale_x_date(date_labels = "%Y-%m-%d") +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = "Date", y = "Frequency", title = "Conformal") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            plot.title = element_text(size = 8))

    p_para / p_conf + plot_annotation(tag_levels = list(c("", "")))
  } else {
    ggplot() + annotate("text", x = .5, y = .5,
                        label = "Interval comparison unavailable") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

# Add panel label manually to the combined plot
fig3c_labeled <- fig3c + plot_annotation(title = NULL) &
  theme(plot.tag = element_blank())
# Wrap so we can use panel_label
fig3c_wrap <- wrap_elements(fig3c) + ggtitle(NULL)

# ── Panel d: Winkler score — Parametric vs Conformal only ────────────────────

fig3d <- tryCatch({
  winkler_data <- map_dfr(calibration$calibration_comparison, function(cc) {
    cc$metrics
  }) |>
    filter(!is.na(winkler), method %in% c("Parametric", "Conformal"))

  ggplot(winkler_data, aes(x = dataset, y = winkler, fill = method)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.5) +
    scale_fill_manual(values = pal_2methods) +
    labs(x = NULL, y = "Winkler score (lower = better)", fill = NULL) +
    panel_label("d") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          legend.position = "bottom",
          legend.key.size = unit(3, "mm"))
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

figure3 <- (fig3a | fig3b) /
            (wrap_elements(fig3c) | fig3d) +
  plot_layout(heights = c(1, 1.6))

save_figure(figure3, "figure3", width_mm = 180, height_mm = 160)

###############################################################################
# FIGURE 4: Decision & Sample Size Impact (180 x 90 mm)
###############################################################################

cat("  [Figure 4] Decision and sample size impact...\n")

# ── Panel a: Vaccine trigger timeline ─────────────────────────────────────────

fig4a <- tryCatch({
  if (nrow(decision) > 0 && !is.null(ba2_data)) {
    # Pick one dataset example (prefer US built-in)
    ba2_lineage <- ba2_data |>
      filter(grepl("BA\\.2|BA2", lineage, ignore.case = TRUE)) |>
      pull(lineage) |> unique()
    if (length(ba2_lineage) == 0) {
      ba2_lineage <- ba2_data |>
        group_by(lineage) |>
        summarise(m = max(proportion, na.rm = TRUE), .groups = "drop") |>
        slice_max(m, n = 1) |> pull(lineage)
    }

    freq_df <- ba2_data |>
      filter(lineage == ba2_lineage[1]) |>
      select(date, frequency = proportion)

    # Get trigger info — Parametric and Conformal only
    triggers <- decision |>
      filter(method %in% c("Parametric", "Conformal"), !is.na(trigger_date)) |>
      slice_head(n = 2)

    actual_date <- triggers$actual_crossing_date[1]

    p <- ggplot(freq_df, aes(x = date, y = frequency)) +
      geom_hline(yintercept = 0.30, linetype = "dashed",
                 linewidth = LW_REF, colour = oi["gray"])

    # Optimal window
    if (!is.na(actual_date)) {
      p <- p + annotate("rect",
                 xmin = actual_date - 7, xmax = actual_date + 7,
                 ymin = -Inf, ymax = Inf,
                 fill = oi["green"], alpha = 0.08)
    }

    p <- p +
      geom_line(linewidth = LW_DATA, colour = "black") +
      geom_point(size = 0.5, colour = "black")

    for (i in seq_len(nrow(triggers))) {
      p <- p + geom_vline(xintercept = triggers$trigger_date[i],
                          colour = pal_2methods[triggers$method[i]],
                          linewidth = LW_DATA * 2)
    }

    p + scale_x_date(date_labels = "%Y-%m-%d") +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = "Date", y = "BA.2 frequency") +
      panel_label("a") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    ggplot() + annotate("text", x = .5, y = .5,
                        label = "Decision data unavailable") + theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel b: Sample size vs MAE and coverage ─────────────────────────────────

fig4b <- tryCatch({
  if (!is.null(sample_size)) {
    ss <- sample_size$summary

    # Long format for faceting
    ss_long <- ss |>
      select(sample_size, MAE = mae_mean, Coverage = cov95_mean) |>
      pivot_longer(cols = c(MAE, Coverage),
                   names_to = "metric", values_to = "value") |>
      mutate(value_pct = value * 100)

    # Reference lines with linetype redundancy for grayscale print
    ref_mae <- tibble(metric = "MAE", yref = 5, lbl = "5%", lt = "dotted")
    ref_cov <- tibble(metric = "Coverage", yref = 90, lbl = "90%", lt = "dashed")
    ref_lines <- bind_rows(ref_mae, ref_cov)

    # Max x for label placement
    x_max <- max(ss_long$sample_size)

    ggplot(ss_long, aes(x = sample_size, y = value_pct)) +
      geom_line(linewidth = LW_DATA, colour = oi["dark_blue"]) +
      geom_point(size = 1.2, colour = oi["dark_blue"]) +
      geom_hline(data = ref_mae, aes(yintercept = yref),
                 linetype = "dotted", linewidth = LW_REF,
                 colour = oi["vermillion"]) +
      geom_hline(data = ref_cov, aes(yintercept = yref),
                 linetype = "dashed", linewidth = LW_REF,
                 colour = oi["vermillion"]) +
      geom_text(data = ref_lines,
                aes(x = x_max, y = yref, label = lbl),
                hjust = 1.1, vjust = -0.5, size = 6 / .pt,
                family = "Arial", colour = oi["vermillion"]) +
      facet_wrap(~metric, scales = "free_y", nrow = 1) +
      scale_x_log10(breaks = c(100, 200, 500, 1000, 2000, 5000)) +
      labs(x = "Sequences per period", y = "%") +
      panel_label("b") +
      theme_nature() +
      theme(strip.text = element_text(size = 8))
  } else {
    ggplot() + annotate("text", x = .5, y = .5,
                        label = "Sample size data unavailable") + theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

figure4 <- fig4a + fig4b + plot_layout(widths = c(1, 1))
save_figure(figure4, "figure4", width_mm = 180, height_mm = 90)

###############################################################################
# FIGURE 5: Advanced Analyses — Supplementary (180 x 140 mm)
###############################################################################

cat("  [Figure 5] Advanced analyses...\n")

# ── Panel a: Fitness decomposition ──────────────────────────────────────────

fig5a <- tryCatch({
  if (!is.null(fitness) && !is.null(fitness$decomposition)) {
    decomp_obj <- fitness$decomposition
    decomp <- if (inherits(decomp_obj, "fitness_decomposition")) {
      decomp_obj$decomposition
    } else if (is.data.frame(decomp_obj)) decomp_obj else NULL
    sens <- fitness$sensitivity

    if (!is.null(decomp) && nrow(decomp) > 0) {
      decomp_long <- decomp |>
        filter(!is.na(transmissibility_fraction)) |>
        select(lineage, Transmissibility = beta,
               Immune_escape = escape_contribution) |>
        pivot_longer(cols = c(Transmissibility, Immune_escape),
                     names_to = "component", values_to = "value")

      p <- ggplot(decomp_long, aes(x = lineage, y = value, fill = component)) +
        geom_col(position = "stack", width = 0.5) +
        scale_fill_manual(values = pal_fitness) +
        labs(x = "Variant", y = expression(delta ~ "(growth advantage)"),
             fill = NULL) +
        panel_label("a") +
        theme_nature() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")

      # Add sensitivity whiskers if available
      if (!is.null(sens) && nrow(sens) > 0) {
        whisker <- sens |>
          filter(abs(perturbation) == 0.20) |>
          group_by(lineage) |>
          summarise(
            esc_lo = min(escape_contribution, na.rm = TRUE),
            esc_hi = max(escape_contribution, na.rm = TRUE),
            .groups = "drop"
          )
        # Total bar height = beta + escape_contribution
        bar_tops <- decomp |>
          filter(!is.na(transmissibility_fraction)) |>
          select(lineage, beta, escape_contribution) |>
          mutate(total = beta + escape_contribution)

        whisker <- whisker |> left_join(bar_tops, by = "lineage")
        if (nrow(whisker) > 0) {
          p <- p + geom_errorbar(
            data = whisker,
            aes(x = lineage,
                ymin = beta + esc_lo,
                ymax = beta + esc_hi),
            inherit.aes = FALSE, width = 0.15, linewidth = LW_DATA
          )
        }
      }
      p
    } else {
      ggplot() + annotate("text", x = .5, y = .5,
                          label = "Fitness data unavailable") + theme_void()
    }
  } else {
    ggplot() + annotate("text", x = .5, y = .5,
                        label = "Fitness data unavailable") + theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel b: Training window vs KS D ────────────────────────────────────────

fig5b <- tryCatch({
  if (!is.null(window_anal) && nrow(window_anal) > 0) {
    trend <- cor(window_anal$window_weeks, window_anal$ks_D,
                 method = "spearman", use = "complete.obs")

    ggplot(window_anal, aes(x = window_weeks, y = ks_D)) +
      geom_line(linewidth = LW_DATA, colour = oi["dark_blue"]) +
      geom_point(size = 1.5, colour = oi["dark_blue"]) +
      annotate("text", x = max(window_anal$window_weeks), y = Inf,
               label = sprintf("Spearman r = %.2f", trend),
               hjust = 1, vjust = 1.5, size = 6 / .pt, family = "Arial") +
      labs(x = "Training window (weeks)", y = "KS D statistic") +
      panel_label("b") +
      theme_nature()
  } else {
    ggplot() + annotate("text", x = .5, y = .5,
                        label = "Window data unavailable") + theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel c: Influenza demonstration ─────────────────────────────────────────

fig5c <- tryCatch({
  flu_nm <- names(influenza)[1]
  flu <- influenza[[flu_nm]]

  if (!is.null(flu$data)) {
    flu_obj <- flu$data

    # lfq_data: .date, .lineage, .freq
    flu_plot <- flu_obj |>
      select(date = .date, subtype = .lineage, frequency = .freq)

    # Top 3 subtypes
    top3 <- flu_plot |>
      group_by(subtype) |>
      summarise(peak = max(frequency, na.rm = TRUE), .groups = "drop") |>
      slice_max(peak, n = 3) |>
      pull(subtype)

    flu_plot <- flu_plot |> filter(subtype %in% top3)

    source_label <- ifelse(flu$source == "simulated",
                           "Influenza (simulated)", "Influenza (WHO FluNet)")

    ggplot(flu_plot, aes(x = date, y = frequency, colour = subtype)) +
      geom_line(linewidth = LW_DATA) +
      scale_x_date(date_labels = "%Y-%m-%d") +
      scale_y_continuous(labels = scales::percent_format()) +
      scale_colour_manual(values = pal_n(length(top3))) +
      labs(x = "Date", y = "Frequency", colour = NULL,
           title = source_label) +
      panel_label("c") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            plot.title = element_text(size = 8))
  } else {
    ggplot() + annotate("text", x = .5, y = .5,
                        label = "No influenza data") + theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel d: EVOI diminishing returns ────────────────────────────────────────

fig5d <- tryCatch({
  if (!is.null(evoi) && nrow(evoi) > 0 && any(!is.na(evoi$evoi))) {
    evoi_plot <- evoi |> filter(!is.na(evoi))

    ggplot(evoi_plot, aes(x = n_additional, y = evoi)) +
      geom_line(linewidth = LW_DATA, colour = oi["dark_blue"]) +
      geom_point(size = 1.5, colour = oi["dark_blue"]) +
      labs(x = "Additional sequences", y = "EVOI (variance reduction)") +
      panel_label("d") +
      theme_nature()
  } else if (!is.null(evoi) && nrow(evoi) > 0 && any(!is.na(evoi$marginal_evoi))) {
    evoi_plot <- evoi |> filter(!is.na(marginal_evoi))

    ggplot(evoi_plot, aes(x = n_additional, y = marginal_evoi)) +
      geom_line(linewidth = LW_DATA, colour = oi["dark_blue"]) +
      geom_point(size = 1.5, colour = oi["dark_blue"]) +
      labs(x = "Additional sequences", y = "Marginal EVOI") +
      panel_label("d") +
      theme_nature()
  } else {
    ggplot() + annotate("text", x = .5, y = .5,
                        label = "EVOI data unavailable") + theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = .5, y = .5, label = paste("Error:", e$message)) +
    theme_void()
})

figure5 <- (fig5a | fig5b) /
            (fig5c | fig5d) +
  plot_layout(heights = c(1, 1.2))

save_figure(figure5, "figure5", width_mm = 180, height_mm = 150)

cat("[08_figures] Complete.\n")
