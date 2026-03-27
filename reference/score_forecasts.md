# Score backtest forecast accuracy

Computes standardized accuracy metrics from backtesting results.

## Usage

``` r
score_forecasts(bt, metrics = c("mae", "rmse", "coverage", "wis"))
```

## Arguments

- bt:

  An `lfq_backtest` object from
  [`backtest()`](https://CuiweiG.github.io/lineagefreq/reference/backtest.md).

- metrics:

  Character vector of metrics to compute:

  - `"mae"`: Mean absolute error of frequency.

  - `"rmse"`: Root mean squared error.

  - `"coverage"`: Proportion within prediction intervals.

  - `"wis"`: Weighted interval score.

## Value

A tibble with columns: `engine`, `horizon`, `metric`, `value`.

## See also

[`compare_models()`](https://CuiweiG.github.io/lineagefreq/reference/compare_models.md)
to rank engines based on scores.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8),
  n_timepoints = 20, seed = 1)
bt <- backtest(sim, engines = "mlr",
  horizons = c(7, 14), min_train = 42)
score_forecasts(bt)
#> # A tibble: 8 × 4
#>   engine horizon metric     value
#>   <chr>    <int> <chr>      <dbl>
#> 1 mlr          7 mae      0.00723
#> 2 mlr          7 rmse     0.00964
#> 3 mlr          7 coverage 0.436  
#> 4 mlr          7 wis      0.00189
#> 5 mlr         14 mae      0.00695
#> 6 mlr         14 rmse     0.00909
#> 7 mlr         14 coverage 0.583  
#> 8 mlr         14 wis      0.00150
# }
```
