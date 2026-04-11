# Extract coefficients from a lineage frequency model

Extract coefficients from a lineage frequency model

## Usage

``` r
# S3 method for class 'lfq_fit'
coef(object, type = c("growth_rate", "all"), ...)
```

## Arguments

- object:

  An `lfq_fit` object.

- type:

  What to return: `"growth_rate"` (default) for growth rates only, or
  `"all"` for intercepts and growth rates.

- ...:

  Ignored.

## Value

Named numeric vector.

## Examples

``` r
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
fit <- fit_model(sim)
coef(fit)
#>          A          B        ref 
#>  0.1793367 -0.2340138  0.0000000 
```
