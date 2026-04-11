# Sequencing power analysis

Estimates the minimum number of sequences needed to detect a lineage at
a given frequency with specified precision.

## Usage

``` r
sequencing_power(target_precision = 0.05, current_freq = 0.02, ci_level = 0.95)
```

## Arguments

- target_precision:

  Desired half-width of the frequency confidence interval. Default 0.05
  (plus/minus 5 percentage points).

- current_freq:

  True or assumed frequency of the target lineage. Can be a vector for
  multiple scenarios. Default 0.02 (2%).

- ci_level:

  Confidence level. Default 0.95.

## Value

A tibble with columns: `current_freq`, `target_precision`, `required_n`,
`ci_level`.

## Details

Uses the normal approximation to the binomial: \$\$n = z^2 \cdot p(1-p)
/ E^2\$\$ where z is the critical value, p is frequency, E is precision.

## Examples

``` r
# How many sequences to estimate a 2% lineage within +/-5%?
sequencing_power()
#> # A tibble: 1 × 4
#>   current_freq target_precision required_n ci_level
#>          <dbl>            <dbl>      <dbl>    <dbl>
#> 1         0.02             0.05         31     0.95

# Multiple scenarios
sequencing_power(current_freq = c(0.01, 0.02, 0.05, 0.10))
#> # A tibble: 4 × 4
#>   current_freq target_precision required_n ci_level
#>          <dbl>            <dbl>      <dbl>    <dbl>
#> 1         0.01             0.05         16     0.95
#> 2         0.02             0.05         31     0.95
#> 3         0.05             0.05         73     0.95
#> 4         0.1              0.05        139     0.95
```
