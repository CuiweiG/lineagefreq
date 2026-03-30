# Print a lineage frequency model

Print a lineage frequency model

## Usage

``` r
# S3 method for class 'lfq_fit'
print(x, ...)
```

## Arguments

- x:

  An `lfq_fit` object.

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## Examples

``` r
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
  n_timepoints = 15, seed = 42)
fit <- fit_model(sim, engine = "mlr")
print(fit)
#> Lineage frequency model (mlr)
#> 3 lineages, 15 time points
#> Date range: 2026-03-30 to 2026-07-06
#> Pivot: "KP.3"
#> 
#> Growth rates (per 7-day unit):
#>   ↑ JN.1: 0.3638
#>   ↑ ref: 0.1122
#> 
#> AIC: 9450; BIC: 9453
```
