# Summarise a lineage frequency model

Summarise a lineage frequency model

## Usage

``` r
# S3 method for class 'lfq_fit'
summary(object, ...)
```

## Arguments

- object:

  An `lfq_fit` object.

- ...:

  Ignored.

## Value

Invisibly returns `object`.

## Examples

``` r
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
  n_timepoints = 15, seed = 42)
fit <- fit_model(sim, engine = "mlr")
summary(fit)
#> Lineage Frequency Model Summary
#> ================================
#> Engine:       mlr 
#> Pivot:        KP.3 
#> Lineages:     3 
#> Time points:  15 
#> Total seqs:   7500 
#> Parameters:   4 
#> Log-lik:      -4721 
#> AIC:          9450 
#> BIC:          9453 
#> 
#> Growth rates (per 7 days):
#> # A tibble: 3 × 6
#>   lineage estimate  lower upper type        pivot
#>   <chr>      <dbl>  <dbl> <dbl> <chr>       <chr>
#> 1 JN.1       0.364 0.338  0.390 growth_rate KP.3 
#> 2 KP.3       0     0      0     growth_rate KP.3 
#> 3 ref        0.112 0.0834 0.141 growth_rate KP.3 
```
