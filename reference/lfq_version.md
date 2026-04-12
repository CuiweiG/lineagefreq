# Package version and system information

Reports lineagefreq version and availability of optional backends.
Useful for reproducibility and bug reports.

## Usage

``` r
lfq_version()
```

## Value

A list with components: `version`, `r_version`, `stan_available`,
`engines`.

## Examples

``` r
lfq_version()
#> $version
#> [1] "0.5.1"
#> 
#> $r_version
#> [1] "4.5.3"
#> 
#> $stan_available
#> [1] FALSE
#> 
#> $engines
#> [1] "mlr"      "hier_mlr" "piantham"
#> 
```
