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

  - `"wis"`: Simplified weighted interval score for the single
    prediction interval stored in the backtest (typically 95%). For the
    full multi-quantile WIS per Bracher et al. (2021), use dedicated
    scoring packages such as 'scoringutils'.

## Value

A tibble with columns: `engine`, `horizon`, `metric`, `value`.

## References

Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
forecasts in an interval format. *PLoS Computational Biology*,
17(2):e1008618.
[doi:10.1371/journal.pcbi.1008618](https://doi.org/10.1371/journal.pcbi.1008618)

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
#> 1 mlr          7 mae      0.00690
#> 2 mlr          7 rmse     0.00932
#> 3 mlr          7 coverage 1      
#> 4 mlr          7 wis      0.00215
#> 5 mlr         14 mae      0.00619
#> 6 mlr         14 rmse     0.00805
#> 7 mlr         14 coverage 1      
#> 8 mlr         14 wis      0.00213
# }
```
