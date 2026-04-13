###############################################################################
# Generate engine comparison figure from saved results
# Run: Rscript inst/extdata/engine_comparison/make_figure.R
###############################################################################

library(lineagefreq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

res <- readRDS("inst/extdata/engine_comparison/engine_comparison_results.rds")

oi <- c(dark_blue = "#0072B2", orange = "#E69F00", green = "#009E73",
        vermillion = "#D55E00", sky_blue = "#56B4E9", gray = "#999999")

eng_colours <- c(MLR = "#0072B2", Piantham = "#E69F00",
                 FGA = "#009E73", GARW = "#D55E00")
eng_order   <- c("MLR", "Piantham", "FGA", "GARW")

LW <- 0.5 * 0.353

theme_fig <- function() {
  theme_classic(base_size = 8, base_family = "Arial") %+replace%
    theme(
      axis.text        = element_text(size = 7, colour = "black"),
      axis.title       = element_text(size = 8),
      strip.text       = element_text(size = 8, face = "bold", hjust = 0),
      strip.background = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line        = element_blank(),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(3, "mm"),
      plot.tag         = element_text(size = 10, face = "bold"),
      plot.margin      = margin(5, 5, 5, 5)
    )
}

# ── Panel a: PIT histograms ─────────────────────────────────────────────

pit_df <- bind_rows(lapply(names(res$calibration), function(eng) {
  cr <- res$calibration[[eng]]
  tibble(engine = toupper(eng), pit = cr$pit_values,
         ks_label = sprintf("KS D = %.3f", cr$ks_D))
})) |>
  mutate(engine = factor(
    ifelse(engine == "PIANTHAM", "Piantham", engine),
    levels = eng_order
  ))

fig_a <- ggplot(pit_df, aes(x = pit)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20,
                 fill = oi["dark_blue"], colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = LW,
             colour = oi["gray"]) +
  geom_text(data = pit_df |> group_by(engine, ks_label) |> slice(1),
            aes(label = ks_label), x = 0.95, y = Inf,
            hjust = 1, vjust = 1.3, size = 5 / .pt, family = "Arial") +
  facet_wrap(~engine, nrow = 1) +
  labs(x = "PIT value", y = "Density", tag = "a") +
  theme_fig()

# ── Panel b: Reliability diagram ────────────────────────────────────────

rel_df <- bind_rows(lapply(names(res$calibration), function(eng) {
  cr <- res$calibration[[eng]]
  cr$reliability |>
    mutate(engine = ifelse(toupper(eng) == "PIANTHAM", "Piantham",
                            toupper(eng)))
})) |>
  mutate(engine = factor(engine, levels = eng_order))

fig_b <- ggplot(rel_df, aes(x = nominal, y = observed,
                             colour = engine, shape = engine)) +
  geom_ribbon(data = tibble(x = seq(0, 1, 0.01)),
              aes(x = x, ymin = pmax(0, x - 0.1), ymax = pmin(1, x + 0.1)),
              inherit.aes = FALSE, fill = "grey90", alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linewidth = LW, colour = "black") +
  geom_line(linewidth = LW * 3) +
  geom_point(size = 1.5) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = eng_colours) +
  scale_shape_manual(values = c(MLR = 16, Piantham = 17, FGA = 15, GARW = 18)) +
  labs(x = "Nominal coverage", y = "Observed coverage",
       colour = NULL, shape = NULL, tag = "b") +
  theme_fig() +
  theme(legend.position = "bottom")

# ── Panel c: KS D bar chart ────────────────────────────────────────────

ks_df <- res$summary |>
  mutate(engine = ifelse(engine == "piantham", "Piantham", toupper(engine)),
         engine = factor(engine, levels = eng_order))

fig_c <- ggplot(ks_df, aes(x = engine, y = ks_D, fill = engine)) +
  geom_col(width = 0.6, colour = NA) +
  scale_fill_manual(values = eng_colours) +
  scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "KS D statistic", tag = "c") +
  theme_fig() +
  theme(legend.position = "none")

# ── Panel d: Coverage at 95% nominal ────────────────────────────────────

cov_df <- bind_rows(lapply(names(res$calibration), function(eng) {
  cr <- res$calibration[[eng]]
  pit <- cr$pit_values
  alpha <- (1 - 0.95) / 2
  obs_cov <- mean(pit >= alpha & pit <= 1 - alpha, na.rm = TRUE)
  tibble(engine = ifelse(eng == "piantham", "Piantham", toupper(eng)),
         coverage = obs_cov)
})) |>
  mutate(engine = factor(engine, levels = eng_order))

fig_d <- ggplot(cov_df, aes(x = engine, y = coverage * 100, fill = engine)) +
  geom_col(width = 0.6, colour = NA) +
  geom_hline(yintercept = 95, linetype = "dashed", linewidth = LW,
             colour = "black") +
  annotate("text", x = 4.4, y = 95, label = "95% nominal",
           hjust = 1, vjust = -0.5, size = 5 / .pt, family = "Arial") +
  scale_fill_manual(values = eng_colours) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Observed 95% coverage (%)", tag = "d") +
  theme_fig() +
  theme(legend.position = "none")

# ── Assemble ────────────────────────────────────────────────────────────

figure <- fig_a / (fig_b | fig_c | fig_d) +
  plot_layout(heights = c(1, 1.2))

dir.create("inst/extdata/engine_comparison", showWarnings = FALSE)
w_in <- 180 / 25.4
h_in <- 140 / 25.4

cairo_pdf("inst/extdata/engine_comparison/figure_engines.pdf",
          width = w_in, height = h_in)
print(figure)
dev.off()

png("inst/extdata/engine_comparison/figure_engines.png",
    width = w_in, height = h_in, units = "in", res = 300, type = "cairo")
print(figure)
dev.off()

cat("Saved figure_engines.pdf and .png\n")
