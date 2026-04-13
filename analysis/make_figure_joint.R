###############################################################################
# make_figure_joint.R — Joint vs marginal conformal prediction visualisation
#
# Produces: submission/figures/fig06_joint_conformal.pdf/.png
# Input:    inst/extdata/engine_comparison/denmark_ba2_collapsed.rds
#
# Run from package root:
#   Rscript analysis/make_figure_joint.R
###############################################################################

Sys.setlocale("LC_TIME", "C")
library(lineagefreq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

setwd("C:/Users/cg223/Desktop/lineagefreq")

Sys.setlocale("LC_TIME", "C")

dir.create("submission/figures", showWarnings = FALSE, recursive = TRUE)

# ── Load Denmark BA.2 collapsed data ─────────────────────────────────────────

dk <- readRDS("inst/extdata/engine_comparison/denmark_ba2_collapsed.rds")
cat("Denmark BA.2:", length(unique(dk$.date)), "dates,",
    length(attr(dk, "lineages")), "lineages\n")

lineages <- attr(dk, "lineages")

# ── Fit model and run both conformal methods ─────────────────────────────────

fit <- fit_model(dk, engine = "mlr")

# Joint conformal
jc <- conformal_forecast_joint(fit, dk, horizon = 14L, ci_level = 0.95,
                                cal_fraction = 0.3, seed = 42)
cat("Joint conformal radius:", round(jc$radius, 3), "\n")
cat("Joint n_cal:", jc$n_cal, "\n")

# Joint calibration from backtest
bt <- readRDS("inst/extdata/engine_comparison/bt_mlr.rds")
jcal <- calibrate_joint(bt)
cat("Joint coverage at 95%:",
    round(jcal$joint_coverage$observed[jcal$joint_coverage$nominal == 0.95] * 100, 1),
    "%\n")

# ── Panel a: Marginal vs joint intervals at a representative origin ──────────
# Choose an origin where BA.2 is between 20-60% for visual clarity

mi <- jc$marginal_intervals
if (nrow(mi) > 0 && ".lower_marginal" %in% names(mi)) {
  # Pick top 3 lineages by range of joint interval width
  top3 <- mi |>
    group_by(.lineage) |>
    summarise(width = mean(.upper_joint - .lower_joint, na.rm = TRUE),
              .groups = "drop") |>
    slice_max(width, n = 3) |>
    pull(.lineage)

  plot_df <- mi |>
    filter(.lineage %in% top3) |>
    pivot_longer(
      cols = c(.lower_joint, .upper_joint, .lower_marginal, .upper_marginal),
      names_to = "bound_type", values_to = "value"
    ) |>
    mutate(
      method = ifelse(grepl("joint", bound_type), "Joint", "Marginal"),
      side   = ifelse(grepl("lower", bound_type), "lower", "upper")
    ) |>
    pivot_wider(names_from = side, values_from = value)

  pal <- c(Marginal = "#0072B2", Joint = "#D55E00")

  fig_a <- ggplot(plot_df, aes(x = .date)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = method),
                alpha = 0.25) +
    geom_point(aes(y = .median), size = 0.8, colour = "black") +
    facet_wrap(~.lineage, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = pal) +
    scale_x_date(date_labels = "%b %Y") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Date", y = "Frequency", fill = NULL, tag = "a") +
    theme_classic(base_size = 8) +
    theme(
      axis.text     = element_text(size = 7, colour = "black"),
      axis.title    = element_text(size = 8),
      strip.text    = element_text(size = 7, face = "bold"),
      panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line     = element_blank(),
      legend.position = "bottom",
      legend.key.size = unit(3, "mm"),
      plot.tag      = element_text(size = 10, face = "bold"),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 5)
    )
} else {
  fig_a <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Marginal comparison not available") +
    theme_void()
}

# ── Panel b: Joint coverage across nominal levels ────────────────────────────

cov_df <- jcal$joint_coverage

fig_b <- ggplot(cov_df, aes(x = nominal, y = observed)) +
  geom_ribbon(data = tibble(x = seq(0, 1, 0.01)),
              aes(x = x, ymin = pmax(0, x - 0.10), ymax = pmin(1, x + 0.10)),
              inherit.aes = FALSE, fill = "grey90", alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25 * 0.353,
              colour = "black") +
  geom_line(colour = "#D55E00", linewidth = 0.5 * 0.353 * 3) +
  geom_point(colour = "#D55E00", size = 1.5) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Nominal coverage", y = "Observed joint coverage", tag = "b") +
  theme_classic(base_size = 8) +
  theme(
    axis.text    = element_text(size = 7, colour = "black"),
    axis.title   = element_text(size = 8),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line    = element_blank(),
    plot.tag     = element_text(size = 10, face = "bold")
  )

# ── Panel c: Multivariate rank histogram ─────────────────────────────────────

rh <- jcal$rank_histogram

fig_c <- ggplot(rh, aes(x = bin, y = density)) +
  geom_col(fill = "#D55E00", colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = rh$expected[1], linetype = "dashed",
             linewidth = 0.25 * 0.353, colour = "grey50") +
  labs(x = "Rank bin", y = "Density", tag = "c") +
  theme_classic(base_size = 8) +
  theme(
    axis.text    = element_text(size = 7, colour = "black"),
    axis.title   = element_text(size = 8),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line    = element_blank(),
    plot.tag     = element_text(size = 10, face = "bold")
  )

# ── Assemble ─────────────────────────────────────────────────────────────────

figure <- fig_a / (fig_b | fig_c) +
  plot_layout(heights = c(1.2, 1))

w <- 180 / 25.4
h <- 130 / 25.4
cairo_pdf("submission/figures/fig06_joint_conformal.pdf", width = w, height = h)
print(figure)
dev.off()

png("submission/figures/fig06_joint_conformal.png",
    width = w, height = h, units = "in", res = 300, type = "cairo")
print(figure)
dev.off()

cat("Saved submission/figures/fig06_joint_conformal.pdf\n")
cat("Saved submission/figures/fig06_joint_conformal.png\n")
