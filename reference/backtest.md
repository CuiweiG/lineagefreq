# Rolling-origin backtesting of lineage frequency models

Evaluates forecast accuracy by repeatedly fitting models on historical
data and comparing predictions to held-out observations. This implements
the evaluation framework described in Abousamra et al. (2024).

## Usage

``` r
backtest(
  data,
  engines = "mlr",
  origins = "weekly",
  horizons = c(7L, 14L, 21L, 28L),
  min_train = 42L,
  ...
)
```

## Arguments

- data:

  An
  [lfq_data](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
  object.

- engines:

  Character vector of engine names to compare. Default `"mlr"`.

- origins:

  How to select forecast origins:

  - `"weekly"` (default): one origin per unique date, starting after
    `min_train` days.

  - An integer: use every Nth date as an origin.

  - A Date vector: use these specific dates as origins.

- horizons:

  Integer vector of forecast horizons in days. Default
  `c(7, 14, 21, 28)`.

- min_train:

  Minimum training window in days. Origins earlier than
  `min(date) + min_train` are skipped. Default 42 (6 weeks).

- ...:

  Additional arguments passed to
  [`fit_model()`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md)
  (e.g., `generation_time` for the Piantham engine).

## Value

An `lfq_backtest` object (tibble subclass) with columns:

- origin_date:

  Date used as the training cutoff.

- target_date:

  Date being predicted.

- horizon:

  Forecast horizon in days.

- engine:

  Engine name.

- lineage:

  Lineage name.

- predicted:

  Predicted frequency (median).

- lower:

  Lower prediction bound.

- upper:

  Upper prediction bound.

- observed:

  Observed frequency at target_date.

## See also

[`score_forecasts()`](https://CuiweiG.github.io/lineagefreq/reference/score_forecasts.md)
to compute accuracy metrics,
[`compare_models()`](https://CuiweiG.github.io/lineagefreq/reference/compare_models.md)
to rank engines.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8),
  n_timepoints = 20, seed = 1)
bt <- backtest(sim, engines = "mlr",
  horizons = c(7, 14), min_train = 42)
bt
#> 
#> ── Backtest results 
#> 75 predictions across 13 origins
#> Engines: "mlr"
#> Horizons: 7, 14 days
#> 
#> # A tibble: 75 × 9
#>    origin_date target_date horizon engine lineage predicted  lower  upper
#>  * <date>      <date>        <int> <chr>  <chr>       <dbl>  <dbl>  <dbl>
#>  1 2026-05-08  2026-05-15        7 mlr    A          0.730  0.699  0.762 
#>  2 2026-05-08  2026-05-15        7 mlr    B          0.0414 0.0324 0.0531
#>  3 2026-05-08  2026-05-15        7 mlr    ref        0.227  0.199  0.259 
#>  4 2026-05-08  2026-05-22       14 mlr    A          0.766  0.730  0.797 
#>  5 2026-05-08  2026-05-22       14 mlr    B          0.0294 0.0222 0.0391
#>  6 2026-05-08  2026-05-22       14 mlr    ref        0.205  0.174  0.240 
#>  7 2026-05-15  2026-05-22        7 mlr    A          0.777  0.748  0.802 
#>  8 2026-05-15  2026-05-22        7 mlr    B          0.0259 0.0196 0.0336
#>  9 2026-05-15  2026-05-22        7 mlr    ref        0.197  0.174  0.224 
#> 10 2026-05-15  2026-05-29       14 mlr    A          0.807  0.779  0.833 
#> # ℹ 65 more rows
#> # ℹ 1 more variable: observed <dbl>
# }
```
