# Compare model engines from backtest scores

Summarises and ranks engines across horizons based on forecast accuracy
scores.

## Usage

``` r
compare_models(scores, by = "engine")
```

## Arguments

- scores:

  Output of
  [`score_forecasts()`](https://CuiweiG.github.io/lineagefreq/reference/score_forecasts.md).

- by:

  Grouping variable(s). Default `"engine"`.

## Value

A tibble with average scores per group, sorted by MAE.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8),
  n_timepoints = 20, seed = 1)
bt <- backtest(sim, engines = "mlr",
  horizons = c(7, 14), min_train = 42)
sc <- score_forecasts(bt)
compare_models(sc)
#> # A tibble: 1 × 5
#>   engine     mae    rmse coverage     wis
#>   <chr>    <dbl>   <dbl>    <dbl>   <dbl>
#> 1 mlr    0.00701 0.00967        1 0.00214
# }
```
