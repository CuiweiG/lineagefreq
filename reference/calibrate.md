# Calibration diagnostics for lineage frequency forecasts

Assesses whether prediction intervals from backtesting are
well-calibrated by computing PIT (Probability Integral Transform)
values, reliability diagrams, and uniformity tests. A perfectly
calibrated forecaster produces PIT values that are uniformly distributed
on \\\[0, 1\]\\.

## Usage

``` r
calibrate(x, observed = NULL, n_bins = 10L)

# S3 method for class 'lfq_backtest'
calibrate(x, observed = NULL, n_bins = 10L)

# S3 method for class 'lfq_forecast'
calibrate(x, observed = NULL, n_bins = 10L)
```

## Arguments

- x:

  An `lfq_backtest` object from
  [`backtest`](https://CuiweiG.github.io/lineagefreq/reference/backtest.md),
  or an `lfq_forecast` object from
  [`forecast`](https://CuiweiG.github.io/lineagefreq/reference/forecast.md)
  together with observed data.

- observed:

  For `lfq_forecast` objects: a numeric vector of observed frequencies
  corresponding to the forecast rows. Ignored for `lfq_backtest` objects
  (which carry their own observed values).

- n_bins:

  Number of bins for the PIT histogram and reliability diagram. Default
  10.

## Value

A `calibration_report` object (S3 class) with components:

- pit_values:

  Numeric vector of PIT values in \\\[0,1\]\\.

- pit_histogram:

  Tibble with `bin`, `count`, `density`, `expected` columns.

- reliability:

  Tibble with `nominal`, `observed` coverage at each level.

- ks_test:

  List with `statistic` and `p_value` from the Kolmogorov-Smirnov test
  for uniformity.

- n:

  Integer; number of forecast-observation pairs used.

## Details

The PIT value for a single forecast-observation pair is defined as the
quantile of the observation within the forecast distribution. Under the
assumption of Gaussian prediction intervals (which the parametric
simulation in
[`forecast`](https://CuiweiG.github.io/lineagefreq/reference/forecast.md)
approximates), the PIT is computed as \\\Phi((y - \hat{y}) /
\hat{\sigma})\\ where \\\hat{y}\\ is the predicted median,
\\\hat{\sigma}\\ is derived from the prediction interval width, and
\\\Phi\\ is the standard normal CDF.

The reliability diagram plots observed coverage against nominal coverage
at levels 10 through 90 percent. Perfect calibration lies on the
diagonal.

## References

Gneiting T, Balabdaoui F, Raftery AE (2007). Probabilistic forecasts,
calibration and sharpness. *Journal of the Royal Statistical Society:
Series B*, 69(2), 243–268.
[doi:10.1111/j.1467-9868.2007.00587.x](https://doi.org/10.1111/j.1467-9868.2007.00587.x)

## See also

[`recalibrate`](https://CuiweiG.github.io/lineagefreq/reference/recalibrate.md)
for post-hoc recalibration,
[`score_forecasts`](https://CuiweiG.github.io/lineagefreq/reference/score_forecasts.md)
for proper scoring rules.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8),
  n_timepoints = 20, seed = 1)
bt <- backtest(sim, engines = "mlr",
  horizons = c(7, 14), min_train = 42)
cal <- calibrate(bt)
#> Warning: ties should not be present for the one-sample Kolmogorov-Smirnov test
cal
#> 
#> ── Calibration report 
#> 75 forecast-observation pairs
#> KS test for PIT uniformity: D = 0.2464, p = 0.000222
#> Mean absolute calibration error: 0.263
# }
```
