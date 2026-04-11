# Tidy an lfq_fit object

Converts model results to a tidy tibble, compatible with the broom
package ecosystem.

## Usage

``` r
tidy.lfq_fit(x, conf.int = TRUE, conf.level = 0.95, ...)
```

## Arguments

- x:

  An `lfq_fit` object.

- conf.int:

  Include confidence intervals? Default `TRUE`.

- conf.level:

  Confidence level. Default 0.95.

- ...:

  Ignored.

## Value

A tibble with columns: `lineage`, `term`, `estimate`, `std.error`,
`conf.low`, `conf.high`.

## Examples

``` r
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
fit <- fit_model(sim)
tidy.lfq_fit(fit)
#> # A tibble: 4 × 6
#>   lineage term        estimate std.error conf.low conf.high
#>   <chr>   <chr>          <dbl>     <dbl>    <dbl>     <dbl>
#> 1 A       intercept     0.0182   0.0491   -0.0779     0.114
#> 2 A       growth_rate   0.179    0.00585   0.168      0.191
#> 3 B       intercept    -0.0318   0.0679   -0.165      0.101
#> 4 B       growth_rate  -0.234    0.0154   -0.264     -0.204
```
