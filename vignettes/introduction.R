## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----setup--------------------------------------------------------------------
library(lineagefreq)

data(sarscov2_us_2022)
head(sarscov2_us_2022)

## ----lfq-data-----------------------------------------------------------------
x <- lfq_data(sarscov2_us_2022,
              lineage = variant,
              date    = date,
              count   = count,
              total   = total)
x

## ----fit----------------------------------------------------------------------
fit <- fit_model(x, engine = "mlr")
fit

## ----growth-advantage---------------------------------------------------------
ga <- growth_advantage(fit,
                       type = "relative_Rt",
                       generation_time = 5)
ga

## ----plot-frequency-----------------------------------------------------------
autoplot(fit, type = "frequency")

## ----plot-advantage-----------------------------------------------------------
autoplot(fit, type = "advantage", generation_time = 5)

## ----forecast-----------------------------------------------------------------
fc <- forecast(fit, horizon = 28)
autoplot(fc)

## ----emergence----------------------------------------------------------------
summarize_emerging(x)

