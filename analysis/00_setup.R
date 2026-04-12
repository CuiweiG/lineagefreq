###############################################################################
# 00_setup.R — System configuration and Nature Methods visual theme
# lineagefreq validation analysis
###############################################################################

cat("[00_setup] Starting...\n")

# ─── System detection ────────────────────────────────────────────────────────

n_cores_available <- parallel::detectCores(logical = TRUE)
n_cores_use       <- max(1L, floor(n_cores_available * 0.80))
ram_gb            <- as.numeric(system("wmic OS get TotalVisibleMemorySize /value", intern = TRUE) |>
                       grep("TotalVisibleMemorySize", x = _, value = TRUE) |>
                       sub("TotalVisibleMemorySize=", "", x = _)) / 1024 / 1024
r_version         <- paste(R.version$major, R.version$minor, sep = ".")
pkg_version       <- as.character(packageVersion("lineagefreq"))

cat(sprintf("  R %s | lineagefreq %s | %d/%d cores (80%%) | %.0f GB RAM\n",
            r_version, pkg_version, n_cores_use, n_cores_available, ram_gb))

# ─── Parallel backend ────────────────────────────────────────────────────────

library(future)
library(furrr)
plan(multisession, workers = n_cores_use)
cat(sprintf("  future::multisession with %d workers\n", n_cores_use))

# ─── Required packages ───────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(lineagefreq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(lubridate)
  library(patchwork)
  library(scales)
  library(kableExtra)
  library(boot)        # for tsboot (block bootstrap)
  library(ISOweek)
})

# ─── Output directories ─────────────────────────────────────────────────────

dir.create("analysis/results", showWarnings = FALSE, recursive = TRUE)
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("analysis/tables",  showWarnings = FALSE, recursive = TRUE)

# ─── Nature Methods visual theme ─────────────────────────────────────────────
# Nature Methods guidelines:
#   - 7 pt axis text, 8 pt axis titles, 9 pt panel/strip titles
#   - Arial font family
#   - Minimal gridlines; theme_classic() base
#   - Line widths: 0.5 pt data lines, 0.25 pt reference lines
#   - Panel border: 0.5 pt solid black

theme_nature <- function(base_size = 8, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Text sizes
      axis.text        = element_text(size = 7, colour = "black"),
      axis.title       = element_text(size = 8, colour = "black"),
      strip.text       = element_text(size = 9, colour = "black", face = "bold",
                                      hjust = 0),
      plot.title       = element_text(size = 9, colour = "black", face = "bold",
                                      hjust = 0, margin = margin(b = 4)),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 8),

      # Panel
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line        = element_blank(),   # border replaces axis lines
      axis.ticks       = element_line(linewidth = 0.25, colour = "black"),

      # Gridlines: none by default
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Strip (facet) background
      strip.background = element_blank(),

      # Legend
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = unit(3, "mm"),
      legend.margin     = margin(0, 0, 0, 0),

      # Spacing
      plot.margin      = margin(5, 5, 5, 5, unit = "pt")
    )
}

# Set as default
theme_set(theme_nature())

# ─── Line width constants (in mm for ggplot2) ────────────────────────────────
# 1 pt = 0.353 mm; ggplot size is in mm
LW_DATA <- 0.5 * 0.353   # 0.5 pt data lines
LW_REF  <- 0.25 * 0.353  # 0.25 pt reference lines

# ─── Unified color palettes ──────────────────────────────────────────────────

pal_countries <- c(
  Denmark     = "#2166AC",
  Germany     = "#4393C3",
  France      = "#92C5DE",
  Netherlands = "#D6604D",
  Spain       = "#B2182B"
)

pal_methods <- c(
  Parametric   = "#2166AC",
  Conformal    = "#B2182B",
  Recalibrated = "#4DAF4A"
)

pal_engines <- c(
  mlr      = "#2166AC",
  piantham = "#D6604D"
)

pal_components <- c(
  Transmissibility = "#2166AC",
  Immune_escape    = "#B2182B"
)

pal_strategies <- c(
  Adaptive      = "#2166AC",
  Proportional  = "#4DAF4A",
  Uniform       = "#636363"
)

# ─── Panel label helper ──────────────────────────────────────────────────────
# Adds Nature Methods-style panel labels (a, b, c, ...) in bold at top-left

add_panel_label <- function(label, x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5) {
  annotate("text", x = x, y = y, label = label,
           fontface = "bold", family = "Arial", size = 9 / .pt,
           hjust = hjust, vjust = vjust)
}

# ─── Colorblind verification ─────────────────────────────────────────────────
# All palettes use RdBu / categorical schemes that are distinguishable under
# deuteranopia and protanopia. Verify with dichromat::dichromat() if available.

verify_colorblind <- function(palette) {
  if (requireNamespace("dichromat", quietly = TRUE)) {
    deutan <- dichromat::dichromat(palette, type = "deutan")
    protan <- dichromat::dichromat(palette, type = "protan")
    # Check that transformed colors remain distinguishable (min pairwise dist > 20)
    check_dist <- function(cols) {
      rgb_mat <- col2rgb(cols)
      d <- as.matrix(dist(t(rgb_mat)))
      diag(d) <- Inf
      min(d)
    }
    d_deutan <- check_dist(deutan)
    d_protan <- check_dist(protan)
    cat(sprintf("  Colorblind check: deutan min_dist=%.0f, protan min_dist=%.0f %s\n",
                d_deutan, d_protan,
                ifelse(min(d_deutan, d_protan) > 15, "[PASS]", "[WARN: low contrast]")))
  } else {
    cat("  dichromat package not available; skipping colorblind verification\n")
  }
}

cat("  Verifying palettes for colorblind safety...\n")
verify_colorblind(pal_countries)
verify_colorblind(pal_methods)

# ─── save_figure() ───────────────────────────────────────────────────────────
# Saves plot as cairo_pdf and 300 dpi PNG
# width_mm / height_mm: figure dimensions in millimetres

save_figure <- function(plot, name, width_mm, height_mm,
                        dir = "analysis/figures") {
  w_in <- width_mm / 25.4
  h_in <- height_mm / 25.4

  pdf_path <- file.path(dir, paste0(name, ".pdf"))
  png_path <- file.path(dir, paste0(name, ".png"))

  cairo_pdf(pdf_path, width = w_in, height = h_in)
  print(plot)
  dev.off()
  cat(sprintf("  Saved: %s (%.0f × %.0f mm)\n", pdf_path, width_mm, height_mm))

  png(png_path, width = w_in, height = h_in, units = "in", res = 300,
      type = "cairo")
  print(plot)
  dev.off()
  cat(sprintf("  Saved: %s (300 dpi)\n", png_path))
}

# ─── save_table() ────────────────────────────────────────────────────────────
# Saves data.frame as LaTeX booktabs table via kableExtra

save_table <- function(df, name, caption = NULL, dir = "analysis/tables") {
  tex_path <- file.path(dir, paste0(name, ".tex"))

  tex <- kableExtra::kbl(df, format = "latex", booktabs = TRUE,
                         caption = caption, linesep = "") |>
    kableExtra::kable_styling(latex_options = c("hold_position"),
                              font_size = 8)

  writeLines(tex, tex_path)
  cat(sprintf("  Saved: %s\n", tex_path))
}

cat("[00_setup] Complete.\n")
