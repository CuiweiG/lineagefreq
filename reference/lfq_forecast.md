# Pipe-friendly forecasting

Pipe-friendly forecasting

## Usage

``` r
lfq_forecast(fit, horizon = 28L, ...)
```

## Arguments

- fit:

  An `lfq_fit` object.

- horizon:

  Forecast horizon in days. Default 28.

- ...:

  Passed to
  [`forecast()`](https://cuiweig.github.io/lineagefreq/reference/forecast.md).

## Value

An `lfq_forecast` object.
