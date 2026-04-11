# Glance at an lfq_fit object

Returns a single-row tibble of model-level summary statistics.

## Usage

``` r
glance.lfq_fit(x, ...)
```

## Arguments

- x:

  An `lfq_fit` object.

- ...:

  Ignored.

## Value

A single-row tibble with columns: `engine`, `n_lineages`,
`n_timepoints`, `nobs`, `df`, `logLik`, `AIC`, `BIC`, `pivot`,
`convergence`.

## Examples

``` r
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
fit <- fit_model(sim)
glance.lfq_fit(fit)
#> # A tibble: 1 × 10
#>   engine n_lineages n_timepoints  nobs    df logLik    AIC    BIC pivot
#>   <chr>       <int>        <int> <int> <int>  <dbl>  <dbl>  <dbl> <chr>
#> 1 mlr             3           20 10000     4 -5645. 11297. 11301. ref  
#> # ℹ 1 more variable: convergence <int>
```
