###############################################################################
# regen_prospective.R — Generate pseudo-prospective evaluation figure
#
# Produces: submission/figures/figure_prospective.pdf
# Input:    cdc_ba2_transition (built-in dataset)
#
# Run from package root:
#   Rscript analysis/regen_prospective.R
###############################################################################

Sys.setlocale("LC_TIME", "C")

devtools::load_all(quiet = TRUE)
library(ggplot2)
library(patchwork)

Sys.setlocale("LC_TIME", "C")

dir.create("submission/figures", showWarnings = FALSE, recursive = TRUE)

dat <- cdc_ba2_transition
lfq <- lfq_data(dat, lineage = lineage, date = date, count = count)
lfq <- collapse_lineages(lfq, min_freq = 0.05)

cat(sprintf("US BA.2: %d dates, %d lineages\n",
            length(unique(lfq$.date)), length(attr(lfq, "lineages"))))

result <- evaluate_prospective(lfq, engine = "mlr", horizons = 14L,
                                min_train = 42L, min_cal = 3L,
                                ci_level = 0.95, gamma = 0.05)
print(result)

p1 <- plot(result, type = "coverage") +
  scale_x_continuous() +
  labs(tag = "a") +
  theme(plot.tag = element_text(size = 10, face = "bold"))

p2 <- plot(result, type = "radius") +
  scale_x_date(date_labels = "%b %Y") +
  labs(tag = "b") +
  theme(plot.tag = element_text(size = 10, face = "bold"))

p3 <- plot(result, type = "comparison") +
  labs(tag = "c") +
  theme(plot.tag = element_text(size = 10, face = "bold"))

fig <- p1 / (p2 | p3)

w <- 180 / 25.4
h <- 120 / 25.4

cairo_pdf("submission/figures/figure_prospective.pdf", width = w, height = h)
print(fig)
dev.off()

cat("Saved submission/figures/figure_prospective.pdf\n")
