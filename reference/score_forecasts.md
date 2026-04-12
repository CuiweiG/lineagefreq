# Score backtest forecast accuracy

Computes standardized accuracy metrics from backtesting results.

## Usage

``` r
score_forecasts(
  bt,
  metrics = c("mae", "rmse", "coverage", "wis", "crps", "log_score", "dss",
    "calibration")
)
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
    prediction interval stored in the backtest (typically 95%).

  - `"crps"`: Continuous Ranked Probability Score, assuming Gaussian
    forecast distribution (Gneiting and Raftery, 2007).

  - `"log_score"`: Logarithmic scoring rule evaluated at the observed
    value under the Gaussian forecast density.

  - `"dss"`: Dawid-Sebastiani Score, a proper scoring rule based on the
    predictive mean and variance.

  - `"calibration"`: Mean squared calibration error across nominal
    coverage levels 10\\

## Value

A tibble with columns: `engine`, `horizon`, `metric`, `value`.

## References

Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
forecasts in an interval format. *PLoS Computational Biology*,
17(2):e1008618.
[doi:10.1371/journal.pcbi.1008618](https://doi.org/10.1371/journal.pcbi.1008618)

Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association*, 102(477), 359–378.
[doi:10.1198/016214506000001437](https://doi.org/10.1198/016214506000001437)

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
#> # A tibble: 16 × 4
#>    engine horizon metric         value
#>    <chr>    <int> <chr>          <dbl>
#>  1 mlr          7 mae          0.00736
#>  2 mlr          7 rmse         0.00985
#>  3 mlr          7 coverage     1      
#>  4 mlr          7 wis          0.00215
#>  5 mlr          7 crps         0.00651
#>  6 mlr          7 log_score   -3.09   
#>  7 mlr          7 dss         -8.01   
#>  8 mlr          7 calibration  0.0744 
#>  9 mlr         14 mae          0.00667
#> 10 mlr         14 rmse         0.00950
#> 11 mlr         14 coverage     1      
#> 12 mlr         14 wis          0.00214
#> 13 mlr         14 crps         0.00639
#> 14 mlr         14 log_score   -3.12   
#> 15 mlr         14 dss         -8.07   
#> 16 mlr         14 calibration  0.0779 
# }
```
