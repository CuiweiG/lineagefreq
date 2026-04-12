#!/usr/bin/env Rscript
# ============================================================
# 00_setup.R — Environment configuration and package loading
# ============================================================
# Run this first. It detects system resources, loads packages,
# sets up the parallel backend, and defines shared utilities
# (themes, palettes, output functions) used by all downstream
# scripts.
# ============================================================

cat("[00_setup] Starting...\n")

# ---- System detection ----
n_cores <- parallel::detectCores()
ram_gb  <- tryCatch({
  raw <- system("powershell -command \"(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB\"",
                intern = TRUE)
  round(as.numeric(raw), 1)
}, error = function(e) NA_real_)

cat("============================================================\n")
cat("System configuration\n")
cat("============================================================\n")
cat(sprintf("  R version:    %s\n", R.version.string))
cat(sprintf("  Platform:     %s\n", R.version$platform))
cat(sprintf("  Cores:        %d\n", n_cores))
cat(sprintf("  RAM:          %.0f GB\n", ram_gb))
cat(sprintf("  Working dir:  %s\n", getwd()))
cat("============================================================\n")

# ---- Package loading ----
# Required for analysis scripts; not part of the lineagefreq package itself.
required <- c("devtools", "ggplot2", "dplyr", "tidyr", "tibble",
              "viridis", "scales", "future", "furrr", "kableExtra",
              "gridExtra")

missing <- required[!vapply(required, requireNamespace,
                            logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  cat(sprintf("Installing missing packages: %s\n",
              paste(missing, collapse = ", ")))
  install.packages(missing)
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(viridis)
library(scales)
library(future)
library(furrr)

# Load the development version of lineagefreq from this repo
devtools::load_all()

pkg_version <- tryCatch(
  as.character(desc::desc_get_version()),
  error = function(e) "unknown"
)
cat(sprintf("  lineagefreq:  %s (dev)\n", pkg_version))

# ---- Parallel backend ----
# R has a hard limit of 128 connections. Each furrr worker uses
# one connection, so we cap at 100 to leave headroom for file I/O.
n_workers <- min(max(floor(n_cores * 0.8), 1L), 100L)
plan(multisession, workers = n_workers)
cat(sprintf("  Workers:      %d (of %d cores)\n", n_workers, n_cores))

# ---- Global ggplot2 theme ----
# Nature-quality standards: clean, minimal, professional.
# Axis text 8pt, titles 10pt, main title 12pt.
# No gridlines except where essential. Classic base.
theme_nature <- function(base_size = 10) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      # Text hierarchy
      axis.text       = element_text(size = 8, colour = "#333333"),
      axis.title      = element_text(size = 10, colour = "#1a1a1a"),
      plot.title       = element_text(size = 12, face = "bold",
                                      colour = "#1a1a1a",
                                      margin = margin(b = 8)),
      plot.subtitle    = element_text(size = 9, colour = "#555555",
                                      margin = margin(b = 6)),
      plot.caption     = element_text(size = 7, colour = "#888888",
                                      hjust = 0),
      strip.text       = element_text(size = 9, face = "bold"),

      # Axes
      axis.line        = element_line(linewidth = 0.4, colour = "#333333"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "#555555"),
      axis.ticks.length = unit(1.5, "mm"),

      # Panel
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid       = element_blank(),

      # Legend
      legend.position  = "bottom",
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(4, "mm"),
      legend.background = element_rect(fill = NA, colour = NA),

      # Margins
      plot.margin = margin(t = 5, r = 10, b = 5, l = 5, unit = "mm")
    )
}
theme_set(theme_nature())

# ---- Color palettes ----
# Engines: muted, distinguishable, colorblind-safe (Wong 2011)
pal_engines <- c(
  mlr      = "#0072B2",   # blue
  hier_mlr = "#E69F00",   # amber
  piantham = "#009E73",   # green
  fga      = "#CC79A7",   # mauve
  garw     = "#56B4E9"    # sky
)

# Lineages: 8-color palette for variant frequency plots
pal_lineages <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733",
  "#999933", "#DDCC77", "#CC6677", "#AA4499"
)

# Sequential palette for heatmaps
pal_heat <- viridis::viridis

# Diverging palette for fitness decomposition
pal_diverge <- c(intrinsic = "#2166AC", escape = "#B2182B")

# ---- Output function ----
# Saves both PDF (for publication) and PNG (for quick viewing).
# Width and height in millimeters, converted to inches for R devices.
save_figure <- function(plot_obj, name,
                        width_mm = 180, height_mm = 100) {
  w_in <- width_mm / 25.4
  h_in <- height_mm / 25.4

  pdf_path <- file.path("analysis", "figures",
                         paste0(name, ".pdf"))
  png_path <- file.path("analysis", "figures",
                         paste0(name, ".png"))

  cairo_pdf(pdf_path, width = w_in, height = h_in)
  print(plot_obj)
  dev.off()

  png(png_path, width = width_mm, height = height_mm,
      units = "mm", res = 300, type = "cairo")
  print(plot_obj)
  dev.off()

  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}

cat("[00_setup] Complete.\n\n")
