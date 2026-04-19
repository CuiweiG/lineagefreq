# Pipe-friendly backtesting + scoring

Runs backtest and returns scores in one step.

## Usage

``` r
lfq_score(
  data,
  engines = "mlr",
  horizons = c(14, 28),
  metrics = c("mae", "coverage"),
  ...
)
```

## Arguments

- data:

  An
  [lfq_data](https://cuiweig.github.io/lineagefreq/reference/lfq_data.md)
  object.

- engines:

  Character vector of engine names.

- horizons:

  Forecast horizons in days.

- metrics:

  Score metrics.

- ...:

  Passed to
  [`backtest()`](https://cuiweig.github.io/lineagefreq/reference/backtest.md).

## Value

A tibble of scores.
