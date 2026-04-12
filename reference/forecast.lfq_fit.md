# Forecast lineage frequencies

Projects lineage frequencies forward in time using the fitted model.
Prediction uncertainty is quantified by parametric simulation from the
estimated parameter distribution.

## Usage

``` r
# S3 method for class 'lfq_fit'
forecast(
  object,
  horizon = 28L,
  ci_level = 0.95,
  n_sim = 1000L,
  sampling_noise = TRUE,
  effective_n = 100L,
  ...
)
```

## Arguments

- object:

  An `lfq_fit` object.

- horizon:

  Number of days to forecast. Default 28 (4 weeks).

- ci_level:

  Confidence level for prediction intervals. Default 0.95.

- n_sim:

  Number of parameter draws for prediction intervals. Default 1000.

- sampling_noise:

  Logical; add multinomial sampling noise to prediction intervals?
  Default `TRUE`. When `TRUE`, each draw includes both parameter
  uncertainty (from the MLE covariance) and observation-level
  multinomial sampling variability.

- effective_n:

  Effective sample size for multinomial sampling noise. Default 100,
  corresponding to a typical weekly sequencing volume per reporting
  unit. Smaller values produce wider (more conservative) prediction
  intervals.

- ...:

  Unused.

## Value

An `lfq_forecast` object (tibble subclass) with columns:

- .date:

  Date.

- .lineage:

  Lineage name.

- .median:

  Median predicted frequency.

- .lower:

  Lower prediction bound.

- .upper:

  Upper prediction bound.

- .type:

  "fitted" or "forecast".

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
fit <- fit_model(sim, engine = "mlr")
fc <- forecast(fit, horizon = 21)
fc
#> 
#> ── Lineage frequency forecast 
#> Engine: mlr
#> Forecast start: 2026-08-24 | Horizon: 21 days
#> CI level: 95%
#> 60 fitted + 9 forecast rows
#> 
#> # A tibble: 69 × 6
#>    .date      .lineage .median .lower .upper .type 
#>  * <date>     <chr>      <dbl>  <dbl>  <dbl> <chr> 
#>  1 2026-04-12 A          0.341  0.341  0.341 fitted
#>  2 2026-04-12 B          0.324  0.324  0.324 fitted
#>  3 2026-04-12 ref        0.335  0.335  0.335 fitted
#>  4 2026-04-19 A          0.408  0.408  0.408 fitted
#>  5 2026-04-19 B          0.257  0.257  0.257 fitted
#>  6 2026-04-19 ref        0.335  0.335  0.335 fitted
#>  7 2026-04-26 A          0.476  0.476  0.476 fitted
#>  8 2026-04-26 B          0.198  0.198  0.198 fitted
#>  9 2026-04-26 ref        0.326  0.326  0.326 fitted
#> 10 2026-05-03 A          0.541  0.541  0.541 fitted
#> # ℹ 59 more rows
# }
```
