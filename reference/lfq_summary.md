# Convert lfq_fit results to a summary tibble

Returns a one-row-per-lineage summary with growth rates, fitted
frequencies at first and last time points, and growth advantage in
multiple scales.

## Usage

``` r
lfq_summary(fit, generation_time = NULL)
```

## Arguments

- fit:

  An `lfq_fit` object.

- generation_time:

  Generation time for Rt calculation.

## Value

A tibble.

## Examples

``` r
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
fit <- fit_model(sim)
lfq_summary(fit, generation_time = 5)
#> # A tibble: 3 × 8
#>   lineage growth_rate gr_lower gr_upper is_pivot relative_Rt Rt_lower Rt_upper
#>   <chr>         <dbl>    <dbl>    <dbl> <lgl>          <dbl>    <dbl>    <dbl>
#> 1 A             0.179    0.168    0.191 FALSE          1.14     1.13     1.15 
#> 2 B            -0.234   -0.264   -0.204 FALSE          0.846    0.828    0.865
#> 3 ref           0        0        0     TRUE           1        1        1    
```
