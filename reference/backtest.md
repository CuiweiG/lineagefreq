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

## Details

Implements the rolling-origin evaluation framework described in
Abousamra et al. (2024), Section 2.4. At each origin date, the model is
fit on data up to that date and forecasts are compared to held-out
future observations. This avoids look-ahead bias and provides an honest
assessment of real-time forecast accuracy.

## References

Abousamra E, Figgins M, Bedford T (2024). Fitness models provide
accurate short-term forecasts of SARS-CoV-2 variant frequency. *PLoS
Computational Biology*, 20(9):e1012443.
[doi:10.1371/journal.pcbi.1012443](https://doi.org/10.1371/journal.pcbi.1012443)

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
#>    origin_date target_date horizon engine lineage predicted lower  upper
#>  * <date>      <date>        <int> <chr>  <chr>       <dbl> <dbl>  <dbl>
#>  1 2026-05-17  2026-05-24        7 mlr    A            0.73  0.64 0.815 
#>  2 2026-05-17  2026-05-24        7 mlr    B            0.04  0.01 0.09  
#>  3 2026-05-17  2026-05-24        7 mlr    ref          0.23  0.14 0.31  
#>  4 2026-05-17  2026-05-31       14 mlr    A            0.77  0.68 0.86  
#>  5 2026-05-17  2026-05-31       14 mlr    B            0.03  0    0.06  
#>  6 2026-05-17  2026-05-31       14 mlr    ref          0.2   0.12 0.3   
#>  7 2026-05-24  2026-05-31        7 mlr    A            0.78  0.68 0.86  
#>  8 2026-05-24  2026-05-31        7 mlr    B            0.02  0    0.0652
#>  9 2026-05-24  2026-05-31        7 mlr    ref          0.19  0.12 0.285 
#> 10 2026-05-24  2026-06-07       14 mlr    A            0.81  0.73 0.88  
#> # ℹ 65 more rows
#> # ℹ 1 more variable: observed <dbl>
# }
```
