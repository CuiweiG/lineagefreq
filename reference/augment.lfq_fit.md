# Augment data with fitted values from an lfq_fit object

Augment data with fitted values from an lfq_fit object

## Usage

``` r
augment.lfq_fit(x, ...)
```

## Arguments

- x:

  An `lfq_fit` object.

- ...:

  Ignored.

## Value

A tibble with columns: `.date`, `.lineage`, `.fitted_freq`, `.observed`,
`.pearson_resid`.

## Examples

``` r
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
fit <- fit_model(sim)
augment.lfq_fit(fit)
#> # A tibble: 60 × 5
#>    .date      .lineage .fitted_freq .observed .pearson_resid
#>    <date>     <chr>           <dbl>     <dbl>          <dbl>
#>  1 2026-03-27 A               0.341     0.332       -0.422  
#>  2 2026-03-27 B               0.324     0.332        0.369  
#>  3 2026-03-27 ref             0.335     0.336        0.0579 
#>  4 2026-04-03 A               0.408     0.442        1.54   
#>  5 2026-04-03 B               0.257     0.248       -0.451  
#>  6 2026-04-03 ref             0.335     0.31        -1.18   
#>  7 2026-04-10 A               0.476     0.472       -0.166  
#>  8 2026-04-10 B               0.198     0.198        0.00250
#>  9 2026-04-10 ref             0.326     0.33         0.175  
#> 10 2026-04-17 A               0.541     0.544        0.137  
#> # ℹ 50 more rows
```
