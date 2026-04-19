# Summarize emerging lineages

Tests whether each lineage's frequency is significantly increasing over
time using a binomial GLM. Useful for early warning of lineages that may
warrant enhanced surveillance.

## Usage

``` r
summarize_emerging(data, threshold = 0.01, min_obs = 3L, p_adjust = "holm")
```

## Arguments

- data:

  An
  [lfq_data](https://cuiweig.github.io/lineagefreq/reference/lfq_data.md)
  object.

- threshold:

  Minimum current frequency to test. Default 0.01.

- min_obs:

  Minimum time points observed. Default 3.

- p_adjust:

  P-value adjustment method. Default `"holm"`.

## Value

A tibble with columns: `lineage`, `first_seen`, `last_seen`,
`n_timepoints`, `current_freq`, `growth_rate`, `p_value`, `p_adjusted`,
`significant`, `direction`.

## Examples

``` r
sim <- simulate_dynamics(
  n_lineages = 4,
  advantages = c(emerging = 1.5, stable = 1.0, declining = 0.7),
  n_timepoints = 12, seed = 42)
summarize_emerging(sim)
#> # A tibble: 2 × 10
#>   lineage  first_seen last_seen  n_timepoints current_freq growth_rate   p_value
#>   <chr>    <date>     <date>            <int>        <dbl>       <dbl>     <dbl>
#> 1 emerging 2026-04-19 2026-07-05           12        0.976      0.0645 0        
#> 2 stable   2026-04-19 2026-07-05           12        0.016     -0.0394 9.90e-101
#> # ℹ 3 more variables: p_adjusted <dbl>, significant <lgl>, direction <chr>
```
