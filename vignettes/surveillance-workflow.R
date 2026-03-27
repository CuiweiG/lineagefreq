## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----load---------------------------------------------------------------------
library(lineagefreq)

data(sarscov2_us_2022)
x <- lfq_data(sarscov2_us_2022,
              lineage = variant,
              date    = date,
              count   = count,
              total   = total)

## ----collapse-----------------------------------------------------------------
x_clean <- collapse_lineages(x, min_freq = 0.02)
attr(x_clean, "lineages")

## ----fit----------------------------------------------------------------------
fit <- fit_model(x_clean, engine = "mlr")
summary(fit)

## ----ga-----------------------------------------------------------------------
ga <- growth_advantage(fit,
                       type = "relative_Rt",
                       generation_time = 5)
ga

## ----ga-plot------------------------------------------------------------------
autoplot(fit, type = "advantage", generation_time = 5)

## ----emerge-------------------------------------------------------------------
emerging <- summarize_emerging(x_clean)
emerging[emerging$significant, ]

## ----forecast-----------------------------------------------------------------
fc <- forecast(fit, horizon = 28)
autoplot(fc)

## ----power--------------------------------------------------------------------
sequencing_power(
  target_precision = 0.05,
  current_freq = c(0.01, 0.02, 0.05, 0.10)
)

## ----tidy---------------------------------------------------------------------
tidy.lfq_fit(fit)
glance.lfq_fit(fit)

