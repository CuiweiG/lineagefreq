###############################################################################
# 08_figures.R — Publication figures: 3 main + 2 extended data
# lineagefreq validation analysis
# Structured per Nature Methods reviewer feedback
###############################################################################

cat("[08_figures] Starting...\n")
source("analysis/00_setup.R")

# English locale for month abbreviations in date axes
tryCatch(Sys.setlocale("LC_TIME", "English"), error = function(e)
  tryCatch(Sys.setlocale("LC_TIME", "en_US.UTF-8"), error = function(e2) NULL))

# ─── Load all results ────────────────────────────────────────────────────────

benchmark   <- readRDS("analysis/results/benchmark_multicountry.rds")
calibration <- readRDS("analysis/results/calibration_comparison.rds")
decision    <- readRDS("analysis/results/decision_impact.rds")
sample_size <- tryCatch(readRDS("analysis/results/sample_size_analysis.rds"),
                        error = function(e) NULL)
window_anal <- tryCatch(readRDS("analysis/results/window_analysis.rds"),
                        error = function(e) NULL)
fitness     <- tryCatch(readRDS("analysis/results/fitness_sensitivity.rds"),
                        error = function(e) NULL)
evoi        <- tryCatch(readRDS("analysis/results/evoi_results.rds"),
                        error = function(e) NULL)
influenza   <- readRDS("analysis/results/influenza_results.rds")

# Built-in data: data.frames with date, lineage, count, proportion
ba2_data <- tryCatch(cdc_ba2_transition, error = function(e) NULL)
jn1_data <- tryCatch(cdc_sarscov2_jn1, error = function(e) NULL)

# ─── Color encoding discipline ──────────────────────────────────────────────
# COLORS   = lineage/variant identity ONLY (Okabe-Ito)
# SHAPES   = method/engine (circle=MLR, triangle=naive, square=conformal)
# GRAYSCALE = dataset identity (dark=US, lighter=European)
# LINETYPE = redundancy for accessibility

oi <- c(orange = "#E69F00", sky_blue = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", dark_blue = "#0072B2", vermillion = "#D55E00",
        pink = "#CC79A7", gray = "#999999")

pal_extended <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9",
                  "#CC79A7", "#F0E442", "#999999", "#000000", "#332288")
pal_n <- function(n) head(pal_extended, max(n, 1))

# Methods: blue vs orange (never red vs green)
pal_methods <- c(Parametric = "#0072B2", Conformal = "#E69F00")

# Fitness components
pal_fitness <- c(Transmissibility = "#0072B2", Immune_escape = "#E69F00")

# Dataset grayscale: US dark, European lighter
pal_ds_gray <- c(US_BA2_builtin = "grey20", US_JN1_builtin = "grey30",
                 Denmark_BA2 = "grey45", France_BA2 = "grey55",
                 Germany_BA2 = "grey65", Netherlands_BA2 = "grey75",
                 Spain_BA2 = "grey85")
shape_ds <- c(US_BA2_builtin = 15, US_JN1_builtin = 15,
              Denmark_BA2 = 16, France_BA2 = 16,
              Germany_BA2 = 16, Netherlands_BA2 = 16, Spain_BA2 = 16)

# ─── Panel label: bold 10pt, one per logical group ──────────────────────────
panel_label <- function(label) {
  annotate("text", x = -Inf, y = Inf, label = label,
           fontface = "bold", family = "Arial", size = 10 / .pt,
           hjust = -0.3, vjust = 1.5)
}

# Placeholder for unavailable data
p_unavail <- function(msg = "Data unavailable") {
  ggplot() + annotate("text", x = .5, y = .5, label = msg) + theme_void()
}

###############################################################################
# MAIN FIGURE 1: Pipeline + Core Data (180 x 130 mm)
###############################################################################

cat("  [Figure 1] Pipeline + core data...\n")

# ── Panel a: Expanded pipeline schematic ─────────────────────────────────────

fig1a <- ggplot() +
  xlim(0, 12) + ylim(0, 6) +

  # Row 1: 4 pipeline boxes at y=4
  # Box 1: Data
  annotate("rect", xmin = 0.1, xmax = 2.4, ymin = 3.5, ymax = 4.5,
           fill = "#E8F0FE", colour = "black", linewidth = 0.3) +
  annotate("text", x = 1.25, y = 4.1, label = "Data",
           family = "Arial", size = 8 / .pt, fontface = "bold") +
  annotate("text", x = 1.25, y = 3.7, label = "lfq_data()",
           family = "Arial", size = 5.5 / .pt, colour = "grey40") +
  # Mini icon: bar chart (raw counts)
  annotate("rect", xmin = 0.7, xmax = 0.9, ymin = 2.8, ymax = 3.2,
           fill = "#0072B2", colour = NA) +
  annotate("rect", xmin = 1.0, xmax = 1.2, ymin = 2.6, ymax = 3.2,
           fill = "#E69F00", colour = NA) +
  annotate("rect", xmin = 1.3, xmax = 1.5, ymin = 2.9, ymax = 3.2,
           fill = "#009E73", colour = NA) +

  # Arrow 1→2
  annotate("segment", x = 2.4, xend = 3.0, y = 4.0, yend = 4.0,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.3) +

  # Box 2: Model
  annotate("rect", xmin = 3.0, xmax = 5.4, ymin = 3.5, ymax = 4.5,
           fill = "#E8F0FE", colour = "black", linewidth = 0.3) +
  annotate("text", x = 4.2, y = 4.1, label = "Model",
           family = "Arial", size = 8 / .pt, fontface = "bold") +
  annotate("text", x = 4.2, y = 3.7, label = "fit_model()",
           family = "Arial", size = 5.5 / .pt, colour = "grey40") +
  # Mini icon: sigmoid curve
  annotate("curve", x = 3.5, y = 2.7, xend = 4.9, yend = 3.2,
           curvature = -0.3, linewidth = 0.4, colour = "#0072B2") +

  # Arrow 2→3
  annotate("segment", x = 5.4, xend = 6.0, y = 4.0, yend = 4.0,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.3) +

  # Box 3: Forecast
  annotate("rect", xmin = 6.0, xmax = 8.4, ymin = 3.5, ymax = 4.5,
           fill = "#E8F0FE", colour = "black", linewidth = 0.3) +
  annotate("text", x = 7.2, y = 4.1, label = "Forecast",
           family = "Arial", size = 8 / .pt, fontface = "bold") +
  annotate("text", x = 7.2, y = 3.7, label = "forecast()",
           family = "Arial", size = 5.5 / .pt, colour = "grey40") +
  # Mini icon: curve with ribbon
  annotate("curve", x = 6.5, y = 2.9, xend = 7.9, yend = 3.15,
           curvature = -0.2, linewidth = 0.4, colour = "#0072B2") +
  annotate("ribbon", x = c(6.5, 7.2, 7.9), ymin = c(2.75, 2.8, 2.9),
           ymax = c(3.05, 3.1, 3.4), stat = "identity",
           fill = "#0072B2", alpha = 0.15) +


  # Arrow 3→4
  annotate("segment", x = 8.4, xend = 9.0, y = 4.0, yend = 4.0,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linewidth = 0.3) +

  # Box 4: Calibrate (dashed = optional)
  annotate("rect", xmin = 9.0, xmax = 11.5, ymin = 3.5, ymax = 4.5,
           fill = "#F5F0FF", colour = "black", linewidth = 0.3,
           linetype = "dashed") +
  annotate("text", x = 10.25, y = 4.1, label = "Calibrate",
           family = "Arial", size = 8 / .pt, fontface = "bold") +
  annotate("text", x = 10.25, y = 3.7, label = "calibrate()",
           family = "Arial", size = 5.5 / .pt, colour = "grey40") +
  # Mini icon: histogram (PIT)
  annotate("rect", xmin = 9.6, xmax = 9.8, ymin = 2.8, ymax = 3.3,
           fill = "#0072B2", colour = NA) +
  annotate("rect", xmin = 9.9, xmax = 10.1, ymin = 2.8, ymax = 3.0,
           fill = "#0072B2", colour = NA) +
  annotate("rect", xmin = 10.2, xmax = 10.4, ymin = 2.8, ymax = 2.95,
           fill = "#0072B2", colour = NA) +
  annotate("rect", xmin = 10.5, xmax = 10.7, ymin = 2.8, ymax = 3.1,
           fill = "#0072B2", colour = NA) +
  annotate("rect", xmin = 10.8, xmax = 11.0, ymin = 2.8, ymax = 3.3,
           fill = "#0072B2", colour = NA) +

  # Pseudocode example
  annotate("text", x = 0.1, y = 1.8,
           label = 'x   <- lfq_data(counts, lineage, date, count)',
           family = "mono", size = 5 / .pt, colour = "grey30", hjust = 0) +
  annotate("text", x = 0.1, y = 1.3,
           label = 'fit <- fit_model(x, engine = "mlr")',
           family = "mono", size = 5 / .pt, colour = "grey30", hjust = 0) +
  annotate("text", x = 0.1, y = 0.8,
           label = "fc  <- forecast(fit, horizon = 28)",
           family = "mono", size = 5 / .pt, colour = "grey30", hjust = 0) +

  panel_label("a") +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

# ── Panel b: BA.2 frequency dynamics ────────────────────────────────────────

fig1b <- tryCatch({
  if (!is.null(ba2_data)) {
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
      scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      labs(x = "Date", y = "Frequency", colour = NULL) +
      panel_label("b") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            legend.key.width = unit(5, "mm"))
  } else p_unavail()
}, error = function(e) p_unavail(e$message))

figure1 <- fig1a + fig1b + plot_layout(widths = c(1, 1))
save_figure(figure1, "figure1", width_mm = 180, height_mm = 130)

###############################################################################
# MAIN FIGURE 2: Performance + Calibration (180 x 180 mm)
###############################################################################

cat("  [Figure 2] Performance + calibration...\n")

# ── Panel a: 5-country small multiples ──────────────────────────────────────

fig2a <- tryCatch({
  ecdc_data <- readRDS("analysis/results/ecdc_prepared.rds")
  ecdc_list <- list()
  for (nm in names(ecdc_data)) {
    if (!grepl("BA2", nm)) next
    obj <- ecdc_data[[nm]]
    country <- sub("_BA2$", "", nm)
    cdf <- obj |>
      select(date = .date, lineage = .lineage, frequency = .freq) |>
      mutate(country = country)
    top4 <- cdf |>
      group_by(lineage) |>
      summarise(peak = max(frequency, na.rm = TRUE), .groups = "drop") |>
      slice_max(peak, n = 4) |> pull(lineage)
    ecdc_list[[nm]] <- cdf |> filter(lineage %in% top4)
  }
  if (length(ecdc_list) > 0) {
    ecdc_df <- bind_rows(ecdc_list)
    n_lin <- n_distinct(ecdc_df$lineage)
    ggplot(ecdc_df, aes(x = date, y = frequency, colour = lineage)) +
      geom_line(linewidth = LW_DATA) +
      facet_wrap(~country, nrow = 1) +
      scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      scale_colour_manual(values = pal_n(n_lin)) +
      labs(x = "Date", y = "Frequency", colour = NULL) +
      panel_label("a") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            strip.text = element_text(size = 7),
            legend.position = "bottom", legend.key.width = unit(4, "mm"))
  } else p_unavail("ECDC data unavailable")
}, error = function(e) p_unavail(e$message))

# ── Panel b: Forest plot with naive baseline ─────────────────────────────────

fig2b <- tryCatch({
  mae_df <- benchmark$metrics |>
    filter(engine == "mlr", horizon == 14) |>
    left_join(benchmark$bootstrap_ci |> filter(engine == "mlr", horizon == 14),
              by = c("dataset", "engine", "horizon")) |>
    mutate(mae_pct = mae * 100, lo_pct = mae_lower * 100,
           hi_pct = mae_upper * 100, model = "MLR (this work)")

  # Naive baseline: last observed carried forward
  naive_mae <- map_dfr(names(benchmark$backtest_results), function(nm) {
    bt <- benchmark$backtest_results[[nm]]
    if (bt$engine != "mlr") return(NULL)
    bt_tbl <- bt$result
    bt_14 <- bt_tbl |> filter(horizon == 14)
    if (nrow(bt_14) == 0) return(NULL)
    min_h <- min(bt_tbl$horizon, na.rm = TRUE)
    origin_obs <- bt_tbl |>
      filter(horizon == min_h) |>
      select(origin_date, lineage, naive_pred = observed)
    joined <- bt_14 |>
      left_join(origin_obs, by = c("origin_date", "lineage")) |>
      filter(!is.na(naive_pred), !is.na(observed))
    if (nrow(joined) == 0) return(NULL)
    tibble(dataset = bt$dataset,
           mae_pct = mean(abs(joined$naive_pred - joined$observed)) * 100)
  })

  plot_df <- mae_df |> select(dataset, mae_pct, lo_pct, hi_pct, model)
  if (nrow(naive_mae) > 0) {
    plot_df <- bind_rows(plot_df,
      naive_mae |> mutate(model = "Naive (persistence)",
                          lo_pct = NA_real_, hi_pct = NA_real_))
  }

  # Highlight US_JN1 row
  jn1_datasets <- plot_df$dataset[grepl("JN1", plot_df$dataset)]

  p <- ggplot(plot_df, aes(x = mae_pct, y = reorder(dataset, -mae_pct),
                           colour = model, shape = model))

  if (length(jn1_datasets) > 0) {
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = match(jn1_datasets[1],
                                   levels(reorder(plot_df$dataset, -plot_df$mae_pct))) - 0.4,
                      ymax = match(jn1_datasets[1],
                                   levels(reorder(plot_df$dataset, -plot_df$mae_pct))) + 0.4,
                      fill = "#FFFFDD", alpha = 0.5)
  }

  p + geom_point(size = 2, position = position_dodge(width = 0.4)) +
    geom_errorbarh(aes(xmin = lo_pct, xmax = hi_pct), height = 0.3,
                   linewidth = LW_DATA, position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = c("MLR (this work)" = "#0072B2",
                                   "Naive (persistence)" = "#999999")) +
    scale_shape_manual(values = c("MLR (this work)" = 16,
                                  "Naive (persistence)" = 2)) +
    labs(x = "MAE at 14-day horizon (%)", y = NULL, colour = NULL, shape = NULL) +
    panel_label("b") +
    theme_nature() +
    theme(legend.position = "bottom", legend.key.size = unit(3, "mm"))
}, error = function(e) p_unavail(e$message))

# ── Panel c: PIT histograms as 6 small multiples ────────────────────────────

fig2c <- tryCatch({
  pit_plot_data <- map_dfr(names(calibration$pit_results), function(nm) {
    pr <- calibration$pit_results[[nm]]
    if (pr$engine != "mlr") return(NULL)
    tibble(dataset = pr$dataset, pit = pr$pit_values,
           label = sprintf("%s\nKS D=%.3f", pr$dataset, pr$ks_D))
  })

  # Limit to 6 panels
  labels_to_show <- unique(pit_plot_data$label)[1:min(6, length(unique(pit_plot_data$label)))]
  pit_plot_data <- pit_plot_data |> filter(label %in% labels_to_show)

  ggplot(pit_plot_data, aes(x = pit)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = "#0072B2", colour = "white", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = LW_REF,
               colour = "#999999") +
    facet_wrap(~label, nrow = 2) +
    labs(x = "PIT value", y = "Density") +
    panel_label("c") +
    theme_nature() +
    theme(strip.text = element_text(size = 6))
}, error = function(e) p_unavail(e$message))

# ── Panel d: Reliability diagram (shape + grayscale for dataset) ─────────────

fig2d <- tryCatch({
  rel_df <- calibration$reliability |>
    filter(engine == "mlr", !is.na(observed_coverage))

  # Assign shape: square for US, circle for European
  rel_df <- rel_df |>
    mutate(
      ds_shape = ifelse(grepl("US|builtin", dataset), "US", "Europe"),
      ds_gray  = case_when(
        grepl("US_BA2|ba2_transition", dataset) ~ "grey20",
        grepl("US_JN1|jn1", dataset) ~ "grey35",
        TRUE ~ "grey60"
      )
    )

  ggplot(rel_df, aes(x = nominal_coverage, y = observed_coverage,
                     group = dataset, shape = ds_shape)) +
    geom_ribbon(data = tibble(x = seq(0, 1, 0.01)),
                aes(x = x, ymin = pmax(0, x - 0.10), ymax = pmin(1, x + 0.10)),
                inherit.aes = FALSE, fill = "grey90", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linewidth = LW_REF, colour = "black") +
    geom_line(aes(colour = dataset), linewidth = LW_DATA) +
    geom_point(aes(colour = dataset), size = 1.2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_shape_manual(values = c(US = 15, Europe = 16)) +
    scale_colour_grey(start = 0.15, end = 0.7) +
    labs(x = "Nominal coverage", y = "Observed coverage",
         colour = NULL, shape = NULL) +
    panel_label("d") +
    theme_nature() +
    theme(legend.position = "bottom", legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 5))
}, error = function(e) p_unavail(e$message))

# ── Panel e: Parametric vs Conformal intervals (2 datasets, 2x2) ────────────

fig2e <- tryCatch({
  cc_keys <- names(calibration$calibration_comparison)
  # Pick a European and US example
  dk_key <- grep("Denmark|dk", cc_keys, ignore.case = TRUE, value = TRUE)
  us_key <- grep("BA2_builtin|ba2_transition", cc_keys, ignore.case = TRUE, value = TRUE)
  if (length(dk_key) == 0) dk_key <- cc_keys[1]
  if (length(us_key) == 0) us_key <- cc_keys[min(2, length(cc_keys))]

  make_interval_subpanel <- function(cc_key, dataset_label, method,
                                     lo_col, hi_col, colour) {
    cal_comp <- calibration$calibration_comparison[[cc_key]]
    fcast_data <- if (method == "Parametric") cal_comp$parametric
                  else cal_comp$conformal
    if (is.null(fcast_data)) return(p_unavail())

    top_lin <- fcast_data |>
      group_by(lineage) |>
      summarise(peak = max(observed, na.rm = TRUE), .groups = "drop") |>
      slice_max(peak, n = 1) |> pull(lineage)
    min_h <- min(fcast_data$horizon, na.rm = TRUE)
    df <- fcast_data |> filter(lineage == top_lin[1], horizon == min_h)
    if (nrow(df) == 0) return(p_unavail())

    outside <- if (lo_col %in% names(df))
      df |> filter(observed < .data[[lo_col]] | observed > .data[[hi_col]])
    else tibble()

    lo <- if (lo_col %in% names(df)) lo_col else "lower"
    hi <- if (hi_col %in% names(df)) hi_col else "upper"

    ggplot(df, aes(x = origin_date)) +
      geom_ribbon(aes(ymin = .data[[lo]], ymax = .data[[hi]]),
                  fill = colour, alpha = 0.2) +
      geom_line(aes(y = predicted), colour = colour, linewidth = LW_DATA) +
      geom_point(aes(y = observed), size = 0.4, colour = "black") +
      {if (nrow(outside) > 0)
        geom_point(data = outside, aes(y = observed),
                   size = 0.8, colour = "#D55E00", shape = 16)} +
      scale_x_date(date_labels = "%b %Y") +
      scale_y_continuous(labels = scales::percent_format(),
                         breaks = c(0, 0.5, 1), limits = c(0, 1)) +
      labs(x = NULL, y = NULL,
           title = paste0(dataset_label, " \u2014 ", method)) +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            plot.title = element_text(size = 7))
  }

  dk_label <- sub("_BA2$", "", sub("__.*", "", dk_key[1]))
  us_label <- sub("_BA2.*|__.*", "", us_key[1])

  p_dk_para <- make_interval_subpanel(dk_key[1], dk_label, "Parametric",
                                       "lower", "upper", pal_methods["Parametric"])
  p_dk_conf <- make_interval_subpanel(dk_key[1], dk_label, "Conformal",
                                       "conf_lower", "conf_upper", pal_methods["Conformal"])
  p_us_para <- make_interval_subpanel(us_key[1], us_label, "Parametric",
                                       "lower", "upper", pal_methods["Parametric"])
  p_us_conf <- make_interval_subpanel(us_key[1], us_label, "Conformal",
                                       "conf_lower", "conf_upper", pal_methods["Conformal"])

  wrap_elements((p_dk_para | p_dk_conf) / (p_us_para | p_us_conf))
}, error = function(e) p_unavail(e$message))

# ── Panel f: Winkler score comparison with paired lines ──────────────────────

fig2f <- tryCatch({
  winkler_data <- map_dfr(calibration$calibration_comparison, function(cc) {
    cc$metrics
  }) |>
    filter(!is.na(winkler), method %in% c("Parametric", "Conformal"))

  # Compute improvement
  winkler_wide <- winkler_data |>
    select(dataset, method, winkler) |>
    pivot_wider(names_from = method, values_from = winkler) |>
    mutate(improvement = (Parametric - Conformal) / Parametric * 100)

  p <- ggplot(winkler_data, aes(x = dataset, y = winkler, fill = method)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.5) +
    scale_fill_manual(values = pal_methods)

  # Paired difference lines
  if (nrow(winkler_wide) > 0) {
    for (i in seq_len(nrow(winkler_wide))) {
      p <- p + annotate("segment",
        x = i - 0.15, xend = i + 0.15,
        y = winkler_wide$Parametric[i], yend = winkler_wide$Conformal[i],
        colour = "grey40", linewidth = LW_REF, linetype = "dotted")
    }
  }

  p + labs(x = NULL, y = "Winkler score (lower = better)", fill = NULL) +
    panel_label("f") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          legend.position = "bottom", legend.key.size = unit(3, "mm"))
}, error = function(e) p_unavail(e$message))

figure2 <- (fig2a | fig2b) /
            (fig2c | fig2d) /
            (wrap_elements(fig2e) | fig2f) +
  plot_layout(heights = c(1, 1, 1.2))

save_figure(figure2, "figure2", width_mm = 180, height_mm = 180)

###############################################################################
# MAIN FIGURE 3: Mechanism + Generalization (180 x 100 mm)
###############################################################################

cat("  [Figure 3] Mechanism + generalization...\n")

# ── Panel a: Fitness decomposition ──────────────────────────────────────────

fig3a <- tryCatch({
  if (!is.null(fitness) && !is.null(fitness$decomposition)) {
    decomp_obj <- fitness$decomposition
    decomp <- if (inherits(decomp_obj, "fitness_decomposition"))
      decomp_obj$decomposition else if (is.data.frame(decomp_obj)) decomp_obj
    sens <- fitness$sensitivity

    if (!is.null(decomp) && nrow(decomp) > 0) {
      # Remove reference lineage (pivot) — not interpretable
      decomp_plot <- decomp |>
        filter(!is.na(transmissibility_fraction),
               !is.na(escape_fraction)) |>
        select(lineage, Transmissibility = beta,
               Immune_escape = escape_contribution) |>
        pivot_longer(cols = c(Transmissibility, Immune_escape),
                     names_to = "component", values_to = "value")

      p <- ggplot(decomp_plot, aes(x = lineage, y = value, fill = component)) +
        geom_col(position = "stack", width = 0.5) +
        scale_fill_manual(values = pal_fitness) +
        labs(x = "Variant (relative to reference)",
             y = expression(delta ~ "(growth advantage)"),
             fill = NULL) +
        panel_label("a") +
        theme_nature() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")

      # Sensitivity whiskers from ±20%
      if (!is.null(sens) && nrow(sens) > 0) {
        whisker <- sens |>
          filter(abs(perturbation) == 0.20) |>
          group_by(lineage) |>
          summarise(esc_lo = min(escape_contribution, na.rm = TRUE),
                    esc_hi = max(escape_contribution, na.rm = TRUE),
                    .groups = "drop")
        bar_tops <- decomp |>
          filter(!is.na(transmissibility_fraction)) |>
          select(lineage, beta)
        whisker <- whisker |> left_join(bar_tops, by = "lineage")
        if (nrow(whisker) > 0) {
          p <- p + geom_errorbar(data = whisker,
            aes(x = lineage, ymin = beta + esc_lo, ymax = beta + esc_hi),
            inherit.aes = FALSE, width = 0.15, linewidth = LW_DATA)
        }
      }

      # BA.2 annotation
      p + annotate("text", x = Inf, y = -Inf,
                    label = "BA.2: consistent with\nLyngse et al. 2022",
                    hjust = 1.05, vjust = -0.3, size = 5 / .pt,
                    family = "Arial", colour = "grey40", fontface = "italic")
    } else p_unavail("Fitness data unavailable")
  } else p_unavail("Fitness data unavailable")
}, error = function(e) p_unavail(e$message))

# ── Panel b: Influenza with forecast performance ─────────────────────────────

fig3b <- tryCatch({
  flu_nm <- names(influenza)[1]
  flu <- influenza[[flu_nm]]

  if (!is.null(flu$data)) {
    flu_obj <- flu$data
    flu_plot <- flu_obj |>
      select(date = .date, subtype = .lineage, frequency = .freq)

    top2 <- flu_plot |>
      group_by(subtype) |>
      summarise(peak = max(frequency, na.rm = TRUE), .groups = "drop") |>
      slice_max(peak, n = 2) |> pull(subtype)
    top3 <- flu_plot |>
      group_by(subtype) |>
      summarise(peak = max(frequency, na.rm = TRUE), .groups = "drop") |>
      slice_max(peak, n = 3) |> pull(subtype)

    flu_plot <- flu_plot |> filter(subtype %in% top3)

    # Linetype redundancy: first=solid, second=dashed, third=dotted
    lt_vals <- setNames(c("solid", "dashed", "dotted")[seq_along(top3)], top3)

    source_label <- ifelse(flu$source == "simulated",
                           "Influenza (simulated)", "Influenza (WHO FluNet)")

    # MAE annotation from backtest
    mae_str <- ""
    if (!is.null(flu$backtest)) {
      flu_mae <- mean(abs(flu$backtest$predicted - flu$backtest$observed),
                       na.rm = TRUE)
      mae_str <- sprintf("14d MAE: %.1f%%", flu_mae * 100)
    }

    p <- ggplot(flu_plot, aes(x = date, y = frequency, colour = subtype,
                              linetype = subtype)) +
      geom_line(linewidth = LW_DATA) +
      scale_colour_manual(values = setNames(pal_n(length(top3)), top3)) +
      scale_linetype_manual(values = lt_vals) +
      scale_x_date(date_labels = "%b %Y") +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = "Date", y = "Frequency", colour = NULL, linetype = NULL,
           title = source_label) +
      panel_label("b") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom", plot.title = element_text(size = 8))

    if (nchar(mae_str) > 0) {
      p <- p + annotate("text", x = Inf, y = Inf, label = mae_str,
                         hjust = 1.05, vjust = 1.5, size = 6 / .pt,
                         family = "Arial", fontface = "bold")
    }
    p
  } else p_unavail("No influenza data")
}, error = function(e) p_unavail(e$message))

figure3 <- fig3a + fig3b + plot_layout(widths = c(1, 1))
save_figure(figure3, "figure3", width_mm = 180, height_mm = 100)

###############################################################################
# EXTENDED DATA FIGURE 1: Decision Impact (180 x 100 mm)
###############################################################################

cat("  [Extended Data Figure 1] Decision impact...\n")

# ── Panel a: Vaccine trigger timeline ────────────────────────────────────────

ext1a <- tryCatch({
  if (nrow(decision) > 0 && !is.null(ba2_data)) {
    ba2_lin <- ba2_data |>
      filter(grepl("BA\\.2|BA2", lineage, ignore.case = TRUE)) |>
      pull(lineage) |> unique()
    if (length(ba2_lin) == 0) {
      ba2_lin <- ba2_data |> group_by(lineage) |>
        summarise(m = max(proportion, na.rm = TRUE), .groups = "drop") |>
        slice_max(m, n = 1) |> pull(lineage)
    }
    freq_df <- ba2_data |>
      filter(lineage == ba2_lin[1]) |>
      select(date, frequency = proportion)

    triggers <- decision |>
      filter(method %in% c("Parametric", "Conformal"), !is.na(trigger_date)) |>
      slice_head(n = 2)
    actual_date <- triggers$actual_crossing_date[1]

    p <- ggplot(freq_df, aes(x = date, y = frequency)) +
      geom_hline(yintercept = 0.30, linetype = "dashed",
                 linewidth = LW_REF, colour = "#999999")

    if (!is.na(actual_date)) {
      p <- p +
        annotate("rect", xmin = actual_date - 7, xmax = actual_date + 7,
                 ymin = -Inf, ymax = Inf, fill = "#56B4E9", alpha = 0.08) +
        annotate("text", x = actual_date, y = Inf,
                 label = "Detection\nwindow", vjust = 1.3,
                 size = 5 / .pt, family = "Arial", colour = "#56B4E9")
    }

    p <- p +
      geom_line(linewidth = LW_DATA, colour = "black") +
      geom_point(size = 0.5, colour = "black")

    for (i in seq_len(nrow(triggers))) {
      p <- p + geom_vline(xintercept = triggers$trigger_date[i],
                          colour = pal_methods[triggers$method[i]],
                          linewidth = LW_DATA * 2)
    }

    p + scale_x_date(date_labels = "%b %Y") +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = "Date", y = "BA.2 frequency") +
      panel_label("a") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else p_unavail("Decision data unavailable")
}, error = function(e) p_unavail(e$message))

# ── Panels b & c: Coverage and MAE vs sample size ────────────────────────────

ext1b <- tryCatch({
  if (!is.null(sample_size)) {
    ss <- sample_size$summary
    x_max <- max(ss$sample_size)

    ggplot(ss, aes(x = sample_size, y = cov95_mean * 100)) +
      geom_line(linewidth = LW_DATA, colour = "#0072B2") +
      geom_point(size = 1.2, colour = "#0072B2") +
      geom_hline(yintercept = 90, linetype = "dashed", linewidth = LW_REF,
                 colour = "#D55E00") +
      annotate("text", x = x_max, y = 90, label = "90%",
               hjust = 1.1, vjust = -0.5, size = 6 / .pt,
               family = "Arial", colour = "#D55E00") +
      scale_x_log10(breaks = c(100, 200, 500, 1000, 2000, 5000)) +
      labs(x = "Sequences per period", y = "95% coverage (%)") +
      panel_label("b") +
      theme_nature()
  } else p_unavail("Sample size data unavailable")
}, error = function(e) p_unavail(e$message))

ext1c <- tryCatch({
  if (!is.null(sample_size)) {
    ss <- sample_size$summary
    x_max <- max(ss$sample_size)

    ggplot(ss, aes(x = sample_size, y = mae_mean * 100)) +
      geom_line(linewidth = LW_DATA, colour = "#0072B2") +
      geom_point(size = 1.2, colour = "#0072B2") +
      geom_hline(yintercept = 5, linetype = "dotted", linewidth = LW_REF,
                 colour = "#D55E00") +
      annotate("text", x = x_max, y = 5, label = "5%",
               hjust = 1.1, vjust = -0.5, size = 6 / .pt,
               family = "Arial", colour = "#D55E00") +
      scale_x_log10(breaks = c(100, 200, 500, 1000, 2000, 5000)) +
      labs(x = "Sequences per period", y = "MAE (%)") +
      panel_label("c") +
      theme_nature()
  } else p_unavail("Sample size data unavailable")
}, error = function(e) p_unavail(e$message))

figure_ext1 <- ext1a | ext1b | ext1c
save_figure(figure_ext1, "figure_ext1", width_mm = 180, height_mm = 100)

###############################################################################
# EXTENDED DATA FIGURE 2: Supplementary Analyses (180 x 100 mm)
###############################################################################

cat("  [Extended Data Figure 2] Supplementary analyses...\n")

# ── Panel a: Training window vs KS D ────────────────────────────────────────

ext2a <- tryCatch({
  if (!is.null(window_anal) && nrow(window_anal) > 0) {
    n_w <- nrow(window_anal)
    trend <- cor(window_anal$window_weeks, window_anal$ks_D,
                 method = "spearman", use = "complete.obs")

    ggplot(window_anal, aes(x = window_weeks, y = ks_D)) +
      geom_line(linewidth = LW_DATA, colour = "#0072B2") +
      geom_point(size = 2, colour = "#0072B2") +
      annotate("text", x = max(window_anal$window_weeks), y = Inf,
               label = sprintf("Spearman r = %.2f, n = %d windows", trend, n_w),
               hjust = 1, vjust = 1.5, size = 6 / .pt, family = "Arial") +
      labs(x = "Training window (weeks)", y = "KS D statistic") +
      panel_label("a") +
      theme_nature()
  } else p_unavail("Window data unavailable")
}, error = function(e) p_unavail(e$message))

# ── Panel b: EVOI diminishing returns ────────────────────────────────────────

ext2b <- tryCatch({
  if (!is.null(evoi) && nrow(evoi) > 0) {
    # Pick best available column
    if ("evoi" %in% names(evoi) && any(!is.na(evoi$evoi))) {
      evoi_plot <- evoi |> filter(!is.na(evoi))
      y_col <- "evoi"
      y_lab <- "Expected variance reduction\nper additional sequence"
    } else if ("marginal_evoi" %in% names(evoi) && any(!is.na(evoi$marginal_evoi))) {
      evoi_plot <- evoi |> filter(!is.na(marginal_evoi))
      y_col <- "marginal_evoi"
      y_lab <- "Marginal EVOI"
    } else {
      return(p_unavail("EVOI data unavailable"))
    }

    # Find diminishing returns threshold (where marginal gain < 10% of max)
    vals <- evoi_plot[[y_col]]
    threshold_idx <- which(vals < max(vals, na.rm = TRUE) * 0.10)[1]

    p <- ggplot(evoi_plot, aes(x = n_additional, y = .data[[y_col]])) +
      geom_line(linewidth = LW_DATA, colour = "#0072B2") +
      geom_point(size = 2, colour = "#0072B2") +
      labs(x = "Additional sequences", y = y_lab) +
      panel_label("b") +
      theme_nature()

    if (!is.na(threshold_idx)) {
      thresh_x <- evoi_plot$n_additional[threshold_idx]
      thresh_y <- vals[threshold_idx]
      p <- p +
        annotate("segment", x = thresh_x, xend = thresh_x + 100,
                 y = thresh_y, yend = thresh_y * 1.5,
                 arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
                 linewidth = 0.3, colour = "#D55E00") +
        annotate("text", x = thresh_x + 110, y = thresh_y * 1.5,
                 label = "Diminishing\nreturns", hjust = 0,
                 size = 5 / .pt, family = "Arial", colour = "#D55E00")
    }
    p
  } else p_unavail("EVOI data unavailable")
}, error = function(e) p_unavail(e$message))

figure_ext2 <- ext2a | ext2b
save_figure(figure_ext2, "figure_ext2", width_mm = 180, height_mm = 100)

cat("[08_figures] Complete.\n")
