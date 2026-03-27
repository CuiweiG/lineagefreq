## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----setup--------------------------------------------------------------------
library(lineagefreq)

## ----mlr----------------------------------------------------------------------
data(sarscov2_us_2022)
x <- lfq_data(sarscov2_us_2022,
              lineage = variant, date = date,
              count = count, total = total)

fit_mlr <- fit_model(x, engine = "mlr")
growth_advantage(fit_mlr, type = "growth_rate")

## ----piantham-----------------------------------------------------------------
fit_pian <- fit_model(x, engine = "piantham",
                      generation_time = 5)
growth_advantage(fit_pian, type = "relative_Rt",
                 generation_time = 5)

## ----glance-------------------------------------------------------------------
dplyr::bind_rows(
  glance.lfq_fit(fit_mlr),
  glance.lfq_fit(fit_pian)
)

## ----backtest-----------------------------------------------------------------
bt <- backtest(x,
  engines = c("mlr", "piantham"),
  horizons = c(7, 14, 21),
  min_train = 56,
  generation_time = 5
)
bt

## ----score--------------------------------------------------------------------
sc <- score_forecasts(bt,
  metrics = c("mae", "coverage"))
sc

## ----compare------------------------------------------------------------------
compare_models(sc, by = c("engine", "horizon"))

## ----plot-backtest------------------------------------------------------------
plot_backtest(sc)

