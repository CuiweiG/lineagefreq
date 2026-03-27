# Extract growth advantage estimates

Computes relative fitness of each lineage from a fitted model. Supports
four output formats for different use cases.

## Usage

``` r
growth_advantage(
  fit,
  type = c("growth_rate", "relative_Rt", "selection_coefficient", "doubling_time"),
  generation_time = NULL,
  ci_level = NULL
)
```

## Arguments

- fit:

  An `lfq_fit` object returned by
  [`fit_model()`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md).

- type:

  Output type:

  - `"growth_rate"` (default): raw growth rate delta per `time_scale`
    days (typically per week).

  - `"relative_Rt"`: relative effective reproduction number. Requires
    `generation_time`.

  - `"selection_coefficient"`: relative Rt minus 1. Requires
    `generation_time`.

  - `"doubling_time"`: days for frequency ratio vs pivot to double
    (positive = growing) or halve (negative = declining).

- generation_time:

  Mean generation time in days. Required for `type = "relative_Rt"` and
  `"selection_coefficient"`. Common values: SARS-CoV-2 approximately 5
  days, influenza approximately 3 days.

- ci_level:

  Confidence level for intervals. Default uses the level from the fitted
  model.

## Value

A tibble with columns:

- lineage:

  Lineage name.

- estimate:

  Point estimate.

- lower:

  Lower confidence bound.

- upper:

  Upper confidence bound.

- type:

  Type of estimate.

- pivot:

  Name of pivot (reference) lineage.

## Examples

``` r
sim <- simulate_dynamics(
  n_lineages = 3,
  advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
  n_timepoints = 15, seed = 42
)
fit <- fit_model(sim, engine = "mlr")

# Growth rates per week
growth_advantage(fit, type = "growth_rate")
#> # A tibble: 3 × 6
#>   lineage estimate  lower upper type        pivot
#>   <chr>      <dbl>  <dbl> <dbl> <chr>       <chr>
#> 1 JN.1       0.364 0.338  0.390 growth_rate KP.3 
#> 2 KP.3       0     0      0     growth_rate KP.3 
#> 3 ref        0.112 0.0834 0.141 growth_rate KP.3 

# Relative Rt (needs generation time)
growth_advantage(fit, type = "relative_Rt", generation_time = 5)
#> # A tibble: 3 × 6
#>   lineage estimate lower upper type        pivot
#>   <chr>      <dbl> <dbl> <dbl> <chr>       <chr>
#> 1 JN.1        1.30  1.27  1.32 relative_Rt KP.3 
#> 2 KP.3        1     1     1    relative_Rt KP.3 
#> 3 ref         1.08  1.06  1.11 relative_Rt KP.3 
```
