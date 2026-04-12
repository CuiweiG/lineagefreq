###############################################################################
# 08_figures.R — Nature Methods 3-figure layout
# lineagefreq validation analysis
###############################################################################

cat("[08_figures] Starting...\n")
source("analysis/00_setup.R")

# Load all results
benchmark   <- readRDS("analysis/results/benchmark_multicountry.rds")
calibration <- readRDS("analysis/results/calibration_comparison.rds")
decision    <- readRDS("analysis/results/decision_impact.rds")
sample_size <- tryCatch(readRDS("analysis/results/sample_size_analysis.rds"), error = function(e) NULL)
window_anal <- tryCatch(readRDS("analysis/results/window_analysis.rds"), error = function(e) NULL)
fitness     <- tryCatch(readRDS("analysis/results/fitness_sensitivity.rds"), error = function(e) NULL)
surveillance <- readRDS("analysis/results/surveillance_simulation.rds")
evoi        <- tryCatch(readRDS("analysis/results/evoi_results.rds"), error = function(e) NULL)
influenza   <- readRDS("analysis/results/influenza_results.rds")

# Built-in data for plotting
# cdc_ba2_transition and cdc_sarscov2_jn1 are data.frames with
# columns: date, lineage, count, proportion
ba2_data <- tryCatch(cdc_ba2_transition, error = function(e) NULL)
jn1_data <- tryCatch(cdc_sarscov2_jn1, error = function(e) NULL)

###############################################################################
# FIGURE 1: Overview + Benchmark (180mm × 220mm)
###############################################################################

cat("  [Figure 1] Overview + Benchmark...\n")

# ── Panel a: Pipeline schematic ──────────────────────────────────────────────

fig1a <- ggplot() +
  xlim(0, 10) + ylim(0, 6) +

  # Boxes for pipeline stages
  # 1. Raw Data
  annotate("rect", xmin = 0.2, xmax = 2.2, ymin = 4.5, ymax = 5.5,
           fill = "#E8F0FE", colour = "black", linewidth = 0.3) +
  annotate("text", x = 1.2, y = 5.0, label = "Raw Counts",
           family = "Arial", size = 7/.pt, fontface = "bold") +
  annotate("text", x = 1.2, y = 4.65, label = "lfq_data()",
           family = "Arial", size = 6/.pt, colour = "grey40") +

  # Arrow 1→2

annotate("segment", x = 2.2, xend = 2.8, y = 5.0, yend = 5.0,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.3) +

  # 2. Model Fit
  annotate("rect", xmin = 2.8, xmax = 4.8, ymin = 4.5, ymax = 5.5,
           fill = "#E8F0FE", colour = "black", linewidth = 0.3) +
  annotate("text", x = 3.8, y = 5.0, label = "MLR Fit",
           family = "Arial", size = 7/.pt, fontface = "bold") +
  annotate("text", x = 3.8, y = 4.65, label = "fit_model()",
           family = "Arial", size = 6/.pt, colour = "grey40") +

  # Arrow 2→3
  annotate("segment", x = 4.8, xend = 5.4, y = 5.0, yend = 5.0,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.3) +

  # 3. Forecast
  annotate("rect", xmin = 5.4, xmax = 7.4, ymin = 4.5, ymax = 5.5,
           fill = "#E8F0FE", colour = "black", linewidth = 0.3) +
  annotate("text", x = 6.4, y = 5.0, label = "Forecast",
           family = "Arial", size = 7/.pt, fontface = "bold") +
  annotate("text", x = 6.4, y = 4.65, label = "forecast()",
           family = "Arial", size = 6/.pt, colour = "grey40") +

  # Arrow 3→4
  annotate("segment", x = 7.4, xend = 8.0, y = 5.0, yend = 5.0,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.3) +

  # 4. Validate
  annotate("rect", xmin = 8.0, xmax = 9.8, ymin = 4.5, ymax = 5.5,
           fill = "#FDE8E8", colour = "black", linewidth = 0.3) +
  annotate("text", x = 8.9, y = 5.0, label = "Validate",
           family = "Arial", size = 7/.pt, fontface = "bold") +
  annotate("text", x = 8.9, y = 4.65, label = "backtest()",
           family = "Arial", size = 6/.pt, colour = "grey40") +

  # Lower row: additional modules
  # Calibrate
  annotate("rect", xmin = 0.5, xmax = 3.0, ymin = 2.5, ymax = 3.5,
           fill = "#E8FEE8", colour = "black", linewidth = 0.3) +
  annotate("text", x = 1.75, y = 3.0, label = "Calibrate",
           family = "Arial", size = 7/.pt, fontface = "bold") +
  annotate("text", x = 1.75, y = 2.65, label = "calibrate() / conformal_forecast()",
           family = "Arial", size = 5/.pt, colour = "grey40") +

  # Decompose
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 2.5, ymax = 3.5,
           fill = "#E8FEE8", colour = "black", linewidth = 0.3) +
  annotate("text", x = 5.0, y = 3.0, label = "Decompose",
           family = "Arial", size = 7/.pt, fontface = "bold") +
  annotate("text", x = 5.0, y = 2.65, label = "fitness_decomposition()",
           family = "Arial", size = 5/.pt, colour = "grey40") +

  # Surveil
  annotate("rect", xmin = 7.0, xmax = 9.5, ymin = 2.5, ymax = 3.5,
           fill = "#E8FEE8", colour = "black", linewidth = 0.3) +
  annotate("text", x = 8.25, y = 3.0, label = "Surveil",
           family = "Arial", size = 7/.pt, fontface = "bold") +
  annotate("text", x = 8.25, y = 2.65, label = "surveillance_value()",
           family = "Arial", size = 5/.pt, colour = "grey40") +

  # Connecting arrows from Validate to lower row
  annotate("segment", x = 8.9, xend = 1.75, y = 4.5, yend = 3.5,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.2, linetype = "dashed", colour = "grey50") +
  annotate("segment", x = 8.9, xend = 5.0, y = 4.5, yend = 3.5,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.2, linetype = "dashed", colour = "grey50") +
  annotate("segment", x = 8.9, xend = 8.25, y = 4.5, yend = 3.5,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.2, linetype = "dashed", colour = "grey50") +

  add_panel_label("a") +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

# ── Panel b: BA.2 frequency dynamics (US CDC data) ──────────────────────────

fig1b <- tryCatch({
  if (!is.null(ba2_data)) {
    # ba2_data is a data.frame with: date, lineage, count, proportion
    df_ba2 <- ba2_data |>
      group_by(date, lineage) |>
      summarise(count = sum(count), .groups = "drop") |>
      group_by(date) |>
      mutate(frequency = count / sum(count)) |>
      ungroup()

    ggplot(df_ba2, aes(x = date, y = frequency, colour = lineage)) +
      geom_point(size = 0.8, alpha = 0.6) +
      geom_line(linewidth = LW_DATA) +
      scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "1 month") +
      scale_y_continuous(labels = percent_format()) +
      labs(x = "Date", y = "Frequency", colour = "Lineage") +
      add_panel_label("b") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            legend.key.width = unit(5, "mm"))
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data unavailable") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel c: 5-country ECDC frequency dynamics ──────────────────────────────

fig1c <- tryCatch({
  ecdc_data <- readRDS("analysis/results/ecdc_prepared.rds")

  # Extract BA.2 period data for plotting
  ecdc_plot_list <- list()
  for (nm in names(ecdc_data)) {
    if (!grepl("BA2", nm)) next
    obj <- ecdc_data[[nm]]
    country <- sub("_BA2$", "", nm)

    # lfq_data objects have columns: .date, .lineage, .freq
    df <- obj |>
      select(date = .date, lineage = .lineage, frequency = .freq) |>
      mutate(country = country)
    ecdc_plot_list[[nm]] <- df
  }

  if (length(ecdc_plot_list) > 0) {
    ecdc_plot_df <- bind_rows(ecdc_plot_list) |>
      filter(frequency > 0.01)  # Remove very rare lineages for clarity

    ggplot(ecdc_plot_df, aes(x = date, y = frequency, colour = lineage)) +
      geom_line(linewidth = LW_DATA) +
      facet_wrap(~country, nrow = 1, scales = "free_x") +
      scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "2 months") +
      scale_y_continuous(labels = percent_format()) +
      labs(x = "Date", y = "Frequency", colour = "Lineage") +
      add_panel_label("c") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            legend.position = "bottom",
            strip.text = element_text(size = 7))
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "ECDC data unavailable") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel d: MAE benchmark heatmap ───────────────────────────────────────────

fig1d <- tryCatch({
  heatmap_data <- benchmark$metrics |>
    filter(!is.na(mae)) |>
    mutate(
      horizon_label = paste0(horizon, "d"),
      mae_pct = mae * 100
    )

  ggplot(heatmap_data, aes(x = horizon_label, y = engine, fill = mae_pct)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", mae_pct)),
              size = 6/.pt, family = "Arial") +
    facet_wrap(~dataset, nrow = 1) +
    scale_fill_gradient(low = "#F7F7F7", high = "#B2182B",
                        name = "MAE (%)") +
    labs(x = "Forecast horizon", y = "Engine") +
    add_panel_label("d") +
    theme_nature() +
    theme(strip.text = element_text(size = 6),
          axis.text.x = element_text(angle = 0))
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel e: MAE by horizon with CI ─────────────────────────────────────────

fig1e <- tryCatch({
  mae_ci <- benchmark$metrics |>
    left_join(benchmark$bootstrap_ci,
              by = c("dataset", "engine", "horizon")) |>
    filter(!is.na(mae), engine == "mlr")

  # Add Bedford Lab reference
  bedford_lines <- benchmark$bedford_ref |>
    filter(!is.na(mae)) |>
    mutate(mae_pct = mae * 100)

  ggplot(mae_ci, aes(x = horizon, y = mae * 100, colour = dataset)) +
    geom_line(linewidth = LW_DATA) +
    geom_point(size = 1.2) +
    geom_ribbon(aes(ymin = mae_lower * 100, ymax = mae_upper * 100, fill = dataset),
                alpha = 0.15, colour = NA) +
    # Bedford Lab reference lines
    geom_hline(data = bedford_lines,
               aes(yintercept = mae_pct),
               linetype = "dashed", linewidth = LW_REF, colour = "grey50") +
    scale_x_continuous(breaks = c(7, 14, 21, 28)) +
    labs(x = "Forecast horizon (days)", y = "MAE (%)",
         colour = "Dataset", fill = "Dataset",
         caption = "Dashed lines: Bedford Lab (Abousamra et al. 2024)") +
    add_panel_label("e") +
    theme_nature() +
    theme(legend.position = "bottom",
          plot.caption = element_text(size = 5, colour = "grey50"))
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel f: Runtime comparison ──────────────────────────────────────────────

fig1f <- tryCatch({
  runtime <- benchmark$runtime |>
    filter(!is.na(elapsed_s))

  ggplot(runtime, aes(x = reorder(paste(dataset, engine, sep = "\n"), elapsed_s),
                      y = elapsed_s, fill = engine)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = pal_engines) +
    labs(x = NULL, y = "Runtime (seconds)", fill = "Engine") +
    add_panel_label("f") +
    theme_nature() +
    theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
          legend.position = "bottom") +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("AMD EPYC 9654, %d cores", n_cores_use),
             hjust = 1.1, vjust = 1.5, size = 5/.pt, colour = "grey50",
             family = "Arial")
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Assemble Figure 1 ───────────────────────────────────────────────────────

figure1 <- (fig1a | fig1b) /
            fig1c /
            (fig1d | fig1e) /
            fig1f +
  plot_layout(heights = c(1.5, 1.2, 1.2, 0.8))

save_figure(figure1, "figure1", width_mm = 180, height_mm = 220)

###############################################################################
# FIGURE 2: Calibration (180mm × 250mm)
###############################################################################

cat("  [Figure 2] Calibration...\n")

# ── Panel a: PIT histograms ─────────────────────────────────────────────────

fig2a <- tryCatch({
  pit_plot_data <- list()
  for (nm in names(calibration$pit_results)) {
    pr <- calibration$pit_results[[nm]]
    pit_plot_data[[nm]] <- tibble(
      dataset = pr$dataset,
      engine  = pr$engine,
      pit     = pr$pit_values,
      label   = paste0(pr$dataset, "\nKS D=",
                        sprintf("%.3f", pr$ks_D),
                        ", p=", sprintf("%.2g", pr$ks_p))
    )
  }
  pit_df <- bind_rows(pit_plot_data) |>
    filter(engine == "mlr")  # Focus on MLR

  # Limit to 6 panels
  datasets_to_show <- unique(pit_df$label)[1:min(6, length(unique(pit_df$label)))]
  pit_df <- pit_df |> filter(label %in% datasets_to_show)

  ggplot(pit_df, aes(x = pit)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = "#2166AC", colour = "white", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = LW_REF,
               colour = "grey40") +
    facet_wrap(~label, nrow = 2, scales = "free_y") +
    labs(x = "PIT value", y = "Density") +
    add_panel_label("a") +
    theme_nature() +
    theme(strip.text = element_text(size = 6))
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel b: Reliability diagram ────────────────────────────────────────────

fig2b <- tryCatch({
  rel_df <- calibration$reliability |>
    filter(!is.na(observed_coverage))

  ggplot(rel_df, aes(x = nominal_coverage, y = observed_coverage,
                     colour = dataset)) +
    # Perfect calibration
    geom_abline(slope = 1, intercept = 0, linetype = "solid",
                linewidth = LW_REF, colour = "black") +
    # ±10% tolerance band
    geom_ribbon(data = tibble(x = seq(0, 1, 0.01)),
                aes(x = x, ymin = x - 0.10, ymax = x + 0.10),
                inherit.aes = FALSE, fill = "grey80", alpha = 0.3) +
    geom_line(linewidth = LW_DATA) +
    geom_point(size = 1) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Nominal coverage", y = "Observed coverage",
         colour = "Dataset") +
    add_panel_label("b") +
    theme_nature() +
    theme(legend.position = "bottom",
          legend.key.width = unit(4, "mm"))
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel c: Three prediction interval types on JN.1 ────────────────────────

fig2c <- tryCatch({
  # Find a JN.1 calibration comparison
  jn1_key <- names(calibration$calibration_comparison) |>
    grep("JN1|jn1|JN\\.1", x = _, value = TRUE)

  if (length(jn1_key) == 0) {
    jn1_key <- names(calibration$calibration_comparison)[1]
  }

  cal_comp <- calibration$calibration_comparison[[jn1_key[1]]]

  # Helper to create panel for one method
  # Data from 03_calibration.R: parametric uses lower/upper,
  # conformal uses conf_lower/conf_upper, recalibrated uses recal_lower/recal_upper
  make_interval_panel <- function(fcast_data, method_name, colour,
                                  lo_col = "lower", hi_col = "upper") {
    if (is.null(fcast_data)) return(NULL)
    fcast_df <- if (is.data.frame(fcast_data)) fcast_data else return(NULL)
    if (!lo_col %in% names(fcast_df) || !hi_col %in% names(fcast_df)) return(NULL)

    # Pick dominant lineage
    top_lineage <- fcast_df |>
      group_by(lineage) |>
      summarise(max_obs = max(observed, na.rm = TRUE), .groups = "drop") |>
      slice_max(max_obs, n = 1) |>
      pull(lineage)

    fcast_lin <- fcast_df |>
      filter(lineage == top_lineage[1], horizon == min(horizon))

    if (nrow(fcast_lin) == 0) return(NULL)

    outside <- fcast_lin |>
      filter(observed < .data[[lo_col]] | observed > .data[[hi_col]])

    ggplot(fcast_lin, aes(x = origin_date)) +
      geom_ribbon(aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]),
                  fill = colour, alpha = 0.2) +
      geom_line(aes(y = predicted), colour = colour, linewidth = LW_DATA) +
      geom_point(aes(y = observed), size = 0.6, colour = "black") +
      {if (nrow(outside) > 0)
        geom_point(data = outside, aes(y = observed),
                   size = 0.8, colour = "red", shape = 16)} +
      scale_x_date(date_labels = "%Y-%m-%d") +
      scale_y_continuous(labels = percent_format()) +
      labs(x = NULL, y = "Frequency", title = method_name) +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            plot.title = element_text(size = 8))
  }

  p_para  <- make_interval_panel(cal_comp$parametric, "Parametric",
               pal_methods["Parametric"], "lower", "upper")
  p_conf  <- make_interval_panel(cal_comp$conformal, "Conformal",
               pal_methods["Conformal"], "conf_lower", "conf_upper")
  p_recal <- make_interval_panel(cal_comp$recalibrated, "Recalibrated",
               pal_methods["Recalibrated"], "recal_lower", "recal_upper")

  panels <- Filter(Negate(is.null), list(p_para, p_conf, p_recal))

  if (length(panels) > 0) {
    wrap_plots(panels, ncol = 1) +
      plot_annotation(tag_levels = list("c"))
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No interval data") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel d: Winkler score comparison ────────────────────────────────────────

fig2d <- tryCatch({
  winkler_data <- map_dfr(calibration$calibration_comparison, function(cc) {
    cc$metrics
  }) |>
    filter(!is.na(winkler))

  ggplot(winkler_data, aes(x = method, y = winkler, fill = method)) +
    geom_col(width = 0.6) +
    facet_wrap(~dataset, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = pal_methods) +
    labs(x = "Method", y = "Winkler Score (lower = better)") +
    add_panel_label("d") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          legend.position = "none",
          strip.text = element_text(size = 6))
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel e: Decision impact — vaccine trigger timeline ──────────────────────

fig2e <- tryCatch({
  if (nrow(decision) > 0 && !is.null(ba2_data)) {
    # ba2_data is a data.frame with: date, lineage, count, proportion
    ba2_lineage <- ba2_data |>
      filter(grepl("BA\\.2|BA2", lineage, ignore.case = TRUE)) |>
      pull(lineage) |>
      unique()
    if (length(ba2_lineage) == 0) {
      ba2_lineage <- ba2_data |>
        group_by(lineage) |> summarise(m = max(proportion, na.rm = TRUE), .groups = "drop") |>
        slice_max(m, n = 1) |> pull(lineage)
    }
    freq_df <- ba2_data |>
      filter(lineage == ba2_lineage[1]) |>
      select(date, frequency = proportion)

    # Trigger dates from decision results
    triggers <- decision |>
      filter(!is.na(trigger_date)) |>
      head(3)  # max 3 methods

    # Optimal trigger window
    actual_date <- unique(triggers$actual_crossing_date)[1]

    p <- ggplot(freq_df, aes(x = date, y = frequency)) +
      # Green optimal zone
      annotate("rect",
               xmin = actual_date - 7, xmax = actual_date + 7,
               ymin = -Inf, ymax = Inf,
               fill = "#4DAF4A", alpha = 0.1) +
      geom_hline(yintercept = 0.30, linetype = "dashed",
                 linewidth = LW_REF, colour = "grey50") +
      geom_line(linewidth = LW_DATA, colour = "black") +
      geom_point(size = 0.6)

    if (nrow(triggers) > 0) {
      for (i in seq_len(nrow(triggers))) {
        p <- p + geom_vline(xintercept = triggers$trigger_date[i],
                            colour = pal_methods[triggers$method[i]],
                            linewidth = LW_DATA, linetype = "solid")
      }
    }

    p + scale_x_date(date_labels = "%Y-%m-%d") +
      scale_y_continuous(labels = percent_format()) +
      labs(x = "Date", y = "BA.2 Frequency",
           caption = "Green zone: optimal trigger window (±7 days)") +
      add_panel_label("e") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.caption = element_text(size = 5, colour = "grey50"))
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Decision data unavailable") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel f: Sample size vs coverage curve ───────────────────────────────────

fig2f <- tryCatch({
  if (!is.null(sample_size)) {
    ss <- sample_size$summary

    # Two y-axes via facet trick: two panels side by side
    ss_long <- ss |>
      select(sample_size, mae_mean, cov95_mean) |>
      pivot_longer(cols = c(mae_mean, cov95_mean),
                   names_to = "metric", values_to = "value") |>
      mutate(
        metric_label = ifelse(metric == "mae_mean", "MAE", "95% Coverage"),
        value_pct = value * 100
      )

    ggplot(ss_long, aes(x = sample_size, y = value_pct)) +
      geom_line(linewidth = LW_DATA, colour = "#2166AC") +
      geom_point(size = 1.2, colour = "#2166AC") +
      facet_wrap(~metric_label, scales = "free_y", nrow = 1) +
      # Reference lines
      geom_hline(data = tibble(metric_label = "MAE", y = 5),
                 aes(yintercept = y), linetype = "dashed",
                 linewidth = LW_REF, colour = "#B2182B") +
      geom_hline(data = tibble(metric_label = "95% Coverage", y = 90),
                 aes(yintercept = y), linetype = "dashed",
                 linewidth = LW_REF, colour = "#B2182B") +
      scale_x_log10(breaks = c(100, 200, 500, 1000, 2000, 5000)) +
      labs(x = "Sequences per period", y = "Percent (%)",
           caption = "MAE saturates at ~500 seq; coverage needs ~2000 seq") +
      add_panel_label("f") +
      theme_nature() +
      theme(plot.caption = element_text(size = 5, colour = "grey50"),
            strip.text = element_text(size = 8))
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Sample size data unavailable") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Assemble Figure 2 ───────────────────────────────────────────────────────

figure2 <- (fig2a | fig2b) /
            (fig2c | fig2d) /
            (fig2e | fig2f) +
  plot_layout(heights = c(1, 1.5, 1))

save_figure(figure2, "figure2", width_mm = 180, height_mm = 250)

###############################################################################
# FIGURE 3: Advanced analyses (180mm × 200mm)
###############################################################################

cat("  [Figure 3] Advanced analyses...\n")

# ── Panel a: Fitness decomposition ──────────────────────────────────────────

fig3a <- tryCatch({
  if (!is.null(fitness) && !is.null(fitness$decomposition)) {
    # fitness$decomposition is a fitness_decomposition S3 object
    # $decomposition tibble has: lineage, observed_advantage, beta,
    #   escape_contribution, transmissibility_fraction, escape_fraction
    decomp_obj <- fitness$decomposition
    decomp <- if (inherits(decomp_obj, "fitness_decomposition")) {
      decomp_obj$decomposition
    } else if (is.data.frame(decomp_obj)) decomp_obj else NULL
    sens <- fitness$sensitivity

    if (!is.null(decomp) && is.data.frame(decomp)) {
      # Pivot for stacked bar using correct column names
      decomp_long <- decomp |>
        filter(!is.na(transmissibility_fraction)) |>
        select(lineage, Transmissibility = beta,
               Immune_escape = escape_contribution) |>
        pivot_longer(cols = c(Transmissibility, Immune_escape),
                     names_to = "component", values_to = "value")

      # Sensitivity whiskers from ±20% perturbation
      whisker_data <- NULL
      if (!is.null(sens) && nrow(sens) > 0) {
        whisker_data <- sens |>
          filter(abs(perturbation) == 0.20) |>
          group_by(lineage) |>
          summarise(
            trans_lo  = min(beta, na.rm = TRUE),
            trans_hi  = max(beta, na.rm = TRUE),
            escape_lo = min(escape_contribution, na.rm = TRUE),
            escape_hi = max(escape_contribution, na.rm = TRUE),
            .groups = "drop"
          )
      }

      p <- ggplot(decomp_long, aes(x = lineage, y = value, fill = component)) +
        geom_col(position = "stack", width = 0.5) +
        scale_fill_manual(values = pal_components,
                          labels = c("Immune escape", "Transmissibility")) +
        labs(x = "Variant", y = "Fitness component (δ)",
             fill = "Component") +
        add_panel_label("a") +
        theme_nature() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")

      p
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Decomposition unavailable") +
        theme_void()
    }
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Fitness data unavailable") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel b: Training window vs PIT uniformity ─────────────────────────────

fig3b <- tryCatch({
  if (!is.null(window_anal) && nrow(window_anal) > 0) {
    ggplot(window_anal, aes(x = window_weeks, y = ks_D)) +
      geom_line(linewidth = LW_DATA, colour = "#2166AC") +
      geom_point(size = 1.5, colour = "#2166AC") +
      labs(x = "Training window (weeks)", y = "KS D statistic",
           caption = "Model misspecification accumulates with window length") +
      add_panel_label("b") +
      theme_nature() +
      theme(plot.caption = element_text(size = 5, colour = "grey50"))
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Window analysis unavailable") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel c: Adaptive surveillance simulation (3 sub-panels) ────────────────

fig3c <- tryCatch({
  sim <- surveillance$raw

  # Left: MAE over time by strategy
  mae_traj <- sim |>
    unnest(mae_trajectory) |>
    group_by(strategy, replicate) |>
    mutate(time_idx = row_number()) |>
    ungroup()

  mae_summary <- mae_traj |>
    group_by(strategy, time_idx) |>
    summarise(
      mae_mean = mean(mae_trajectory, na.rm = TRUE),
      mae_lo   = quantile(mae_trajectory, 0.025, na.rm = TRUE),
      mae_hi   = quantile(mae_trajectory, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  p_mae <- ggplot(mae_summary, aes(x = time_idx, y = mae_mean * 100,
                                    colour = strategy, fill = strategy)) +
    geom_ribbon(aes(ymin = mae_lo * 100, ymax = mae_hi * 100), alpha = 0.15,
                colour = NA) +
    geom_line(linewidth = LW_DATA) +
    scale_colour_manual(values = pal_strategies) +
    scale_fill_manual(values = pal_strategies) +
    labs(x = "Time index", y = "MAE (%)") +
    theme_nature() +
    theme(legend.position = "none")

  # Center: Detection delay boxplot
  delay_df <- sim |>
    unnest(detection_delay) |>
    filter(!is.na(detection_delay))

  p_delay <- ggplot(delay_df, aes(x = strategy, y = detection_delay,
                                   fill = strategy)) +
    geom_boxplot(width = 0.5, outlier.size = 0.5) +
    scale_fill_manual(values = pal_strategies) +
    labs(x = NULL, y = "Detection delay (days)") +
    theme_nature() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  # Right: EVOI curve
  p_evoi <- if (!is.null(evoi) && nrow(evoi) > 0 && any(!is.na(evoi$evoi))) {
    ggplot(evoi |> filter(!is.na(evoi)), aes(x = n_sequences, y = evoi)) +
      geom_line(linewidth = LW_DATA, colour = "#2166AC") +
      geom_point(size = 1.2, colour = "#2166AC") +
      labs(x = "Additional sequences", y = "EVOI") +
      theme_nature()
  } else {
    ggplot(evoi |> filter(!is.na(mae_expected)),
           aes(x = n_sequences, y = mae_expected * 100)) +
      geom_line(linewidth = LW_DATA, colour = "#2166AC") +
      geom_point(size = 1.2, colour = "#2166AC") +
      labs(x = "Additional sequences", y = "Expected MAE (%)") +
      theme_nature()
  }

  (p_mae | p_delay | p_evoi) +
    plot_annotation(tag_levels = list("c"))
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel d: Influenza demonstration ─────────────────────────────────────────

fig3d <- tryCatch({
  # Pick first influenza dataset
  flu_nm <- names(influenza)[1]
  flu <- influenza[[flu_nm]]

  if (!is.null(flu$data)) {
    flu_obj <- flu$data
    # flu_obj is lfq_data (tibble with .date, .lineage, .freq)
    flu_plot <- flu_obj |>
      select(date = .date, subtype = .lineage, frequency = .freq)

    source_label <- ifelse(flu$source == "simulated",
                           "SIMULATED DATA", "WHO FluNet")

    ggplot(flu_plot, aes(x = date, y = frequency, colour = subtype)) +
      geom_line(linewidth = LW_DATA) +
      scale_x_date(date_labels = "%Y-%m-%d") +
      scale_y_continuous(labels = percent_format()) +
      labs(x = "Date", y = "Frequency", colour = "Subtype",
           title = source_label) +
      add_panel_label("d") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No influenza data") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Panel e: Identifiability sensitivity ────────────────────────────────────

fig3e <- tryCatch({
  if (!is.null(fitness) && !is.null(fitness$sensitivity) &&
      nrow(fitness$sensitivity) > 0) {
    sens <- fitness$sensitivity

    ggplot(sens, aes(x = perturbation * 100, y = escape_fraction,
                     colour = lineage)) +
      geom_line(linewidth = LW_DATA) +
      geom_point(size = 1) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = LW_REF) +
      labs(x = "Immunity perturbation (%)",
           y = "Immune escape fraction",
           colour = "Lineage") +
      add_panel_label("e") +
      theme_nature() +
      theme(legend.position = "bottom")
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5,
                        label = "Sensitivity data unavailable") +
      theme_void()
  }
}, error = function(e) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
    theme_void()
})

# ── Assemble Figure 3 ───────────────────────────────────────────────────────

figure3 <- (fig3a | fig3b) /
            fig3c /
            (fig3d | fig3e) +
  plot_layout(heights = c(1, 1.2, 1))

save_figure(figure3, "figure3", width_mm = 180, height_mm = 200)

cat("[08_figures] Complete.\n")
