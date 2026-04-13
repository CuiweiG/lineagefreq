###############################################################################
# make_figure_aci_sensitivity.R — ACI learning rate sensitivity analysis
#
# Produces: submission/figures/figS_aci_sensitivity.pdf/.png
# Input:    cdc_ba2_transition (built-in dataset)
#
# Run from package root:
#   Rscript analysis/make_figure_aci_sensitivity.R
# Or in RStudio:
#   source("analysis/make_figure_aci_sensitivity.R")
#
# Expected runtime: ~5-10 minutes (5 gamma values, sequential)
# Requires: lineagefreq >= 0.6.0 installed locally
#   Install with: devtools::install("C:/Users/cg223/Desktop/lineagefreq")
###############################################################################

Sys.setlocale("LC_TIME", "C")
library(lineagefreq)
library(ggplot2)
library(dplyr)

Sys.setlocale("LC_TIME", "C")

pkg_root <- "C:/Users/cg223/Desktop/lineagefreq"
setwd(pkg_root)
dir.create("submission/figures", showWarnings = FALSE, recursive = TRUE)

cat("lineagefreq version:", as.character(packageVersion("lineagefreq")), "\n")

# ── Load and prepare data ────────────────────────────────────────────────────

dat <- cdc_ba2_transition
lfq <- lfq_data(dat, lineage = lineage, date = date, count = count)
lfq <- collapse_lineages(lfq, min_freq = 0.05)
cat(sprintf("US BA.2: %d dates, %d lineages\n",
            length(unique(lfq$.date)), length(attr(lfq, "lineages"))))

# ── Run evaluate_prospective with different gamma values (sequential) ────────

gammas <- c(0.01, 0.02, 0.05, 0.10, 0.20)
results <- list()

for (i in seq_along(gammas)) {
  g <- gammas[i]
  cat(sprintf("[%d/%d] gamma = %.2f ...\n", i, length(gammas), g))

  ep <- evaluate_prospective(
    lfq,
    engine    = "mlr",
    horizons  = 14L,
    min_train = 42L,
    min_cal   = 3L,
    ci_level  = 0.95,
    gamma     = g
  )

  results[[i]] <- tibble(
    gamma           = g,
    coverage_param  = ep$summary$coverage[ep$summary$method == "Parametric"],
    coverage_static = ep$summary$coverage[ep$summary$method == "Static conformal"],
    coverage_aci    = ep$summary$coverage[ep$summary$method == "ACI"],
    width_static    = ep$summary$mean_width[ep$summary$method == "Static conformal"],
    width_aci       = ep$summary$mean_width[ep$summary$method == "ACI"],
    winkler_static  = ep$summary$winkler_score[ep$summary$method == "Static conformal"],
    winkler_aci     = ep$summary$winkler_score[ep$summary$method == "ACI"]
  )
}

results_df <- bind_rows(results)
cat("\nResults:\n")
print(results_df)

# ── Panel a: Coverage vs gamma ───────────────────────────────────────────────

pal <- c("Static conformal" = "#0072B2", "ACI" = "#009E73")

cov_long <- results_df |>
  select(gamma, `Static conformal` = coverage_static,
         ACI = coverage_aci) |>
  tidyr::pivot_longer(cols = c(`Static conformal`, ACI),
                      names_to = "method", values_to = "coverage")

fig_a <- ggplot(cov_long, aes(x = gamma, y = coverage, colour = method,
                               shape = method)) +
  geom_line(linewidth = 0.5 * 0.353 * 3) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.25 * 0.353,
             colour = "grey40") +
  annotate("text", x = max(gammas), y = 0.95, label = "95% target",
           hjust = 1, vjust = -0.5, size = 5 / .pt, colour = "grey40") +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = c("Static conformal" = 16, ACI = 17)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0.85, 1.0)) +
  scale_x_continuous(breaks = gammas) +
  labs(x = expression(paste("ACI learning rate ", gamma)),
       y = "Coverage at 95% nominal",
       colour = NULL, shape = NULL, tag = "a") +
  theme_classic(base_size = 8) +
  theme(
    axis.text    = element_text(size = 7, colour = "black"),
    axis.title   = element_text(size = 8),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line    = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    plot.tag     = element_text(size = 10, face = "bold")
  )

# ── Panel b: Mean width vs gamma ────────────────────────────────────────────

width_long <- results_df |>
  select(gamma, `Static conformal` = width_static,
         ACI = width_aci) |>
  tidyr::pivot_longer(cols = c(`Static conformal`, ACI),
                      names_to = "method", values_to = "width")

fig_b <- ggplot(width_long, aes(x = gamma, y = width, colour = method,
                                 shape = method)) +
  geom_line(linewidth = 0.5 * 0.353 * 3) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = c("Static conformal" = 16, ACI = 17)) +
  scale_x_continuous(breaks = gammas) +
  labs(x = expression(paste("ACI learning rate ", gamma)),
       y = "Mean interval width",
       colour = NULL, shape = NULL, tag = "b") +
  theme_classic(base_size = 8) +
  theme(
    axis.text    = element_text(size = 7, colour = "black"),
    axis.title   = element_text(size = 8),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line    = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    plot.tag     = element_text(size = 10, face = "bold")
  )

# ── Assemble and save ────────────────────────────────────────────────────────

figure <- fig_a | fig_b

w <- 180 / 25.4
h <- 80 / 25.4

cairo_pdf("submission/figures/figS_aci_sensitivity.pdf", width = w, height = h)
print(figure)
dev.off()

png("submission/figures/figS_aci_sensitivity.png",
    width = w, height = h, units = "in", res = 300, type = "cairo")
print(figure)
dev.off()

# Save results for reference
dir.create("inst/extdata/engine_comparison", showWarnings = FALSE, recursive = TRUE)
saveRDS(results_df, "inst/extdata/engine_comparison/aci_sensitivity_results.rds")

cat("Saved submission/figures/figS_aci_sensitivity.pdf\n")
cat("Saved submission/figures/figS_aci_sensitivity.png\n")
cat("Saved inst/extdata/engine_comparison/aci_sensitivity_results.rds\n")
