# Simulate lineage frequency dynamics

Generates synthetic lineage frequency data under a multinomial sampling
model with configurable growth advantages. Useful for model validation,
power analysis, and teaching.

## Usage

``` r
simulate_dynamics(
  n_lineages = 4L,
  n_timepoints = 20L,
  total_per_tp = 500L,
  advantages = NULL,
  reference_name = "ref",
  start_date = Sys.Date(),
  interval = 7L,
  overdispersion = NULL,
  seed = NULL
)
```

## Arguments

- n_lineages:

  Number of lineages including reference. Default 4.

- n_timepoints:

  Number of time points. Default 20.

- total_per_tp:

  Sequences per time point. A single integer (constant across time) or a
  vector of length `n_timepoints` (variable sampling effort). Default
  500.

- advantages:

  Named numeric vector of per-week multiplicative growth advantages for
  non-reference lineages. Length must equal `n_lineages - 1`. Values \>
  1 mean growing faster than reference, \< 1 means declining. If `NULL`,
  random values in (0.8, 1.5).

- reference_name:

  Name of the reference lineage. Default `"ref"`.

- start_date:

  Start date for the time series. Default today.

- interval:

  Days between consecutive time points. Default 7 (weekly data).

- overdispersion:

  If `NULL` (default), standard multinomial sampling. If a positive
  number, Dirichlet-Multinomial sampling with concentration = 1 /
  overdispersion. Larger values = more overdispersion.

- seed:

  Random seed for reproducibility.

## Value

An
[lfq_data](https://cuiweig.github.io/lineagefreq/reference/lfq_data.md)
object with an additional `true_freq` column containing the true
(pre-sampling) frequencies.

## Examples

``` r
# JN.1 grows 1.3x per week, KP.3 declines at 0.9x
sim <- simulate_dynamics(
  n_lineages = 3,
  advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
  n_timepoints = 15,
  seed = 42
)
sim
#> 
#> ── Lineage frequency data 
#> 3 lineages, 15 time points
#> Date range: 2026-04-19 to 2026-07-26
#> Lineages: "JN.1, KP.3, ref"
#> 
#> # A tibble: 45 × 7
#>    .date      .lineage .count true_freq .total .freq .reliable
#>  * <date>     <chr>     <int>     <dbl>  <int> <dbl> <lgl>    
#>  1 2026-04-19 JN.1        156     0.333    500 0.312 TRUE     
#>  2 2026-04-19 KP.3        184     0.333    500 0.368 TRUE     
#>  3 2026-04-19 ref         160     0.333    500 0.32  TRUE     
#>  4 2026-04-26 JN.1        204     0.406    500 0.408 TRUE     
#>  5 2026-04-26 KP.3        143     0.281    500 0.286 TRUE     
#>  6 2026-04-26 ref         153     0.312    500 0.306 TRUE     
#>  7 2026-05-03 JN.1        240     0.483    500 0.48  TRUE     
#>  8 2026-05-03 KP.3        120     0.231    500 0.24  TRUE     
#>  9 2026-05-03 ref         140     0.286    500 0.28  TRUE     
#> 10 2026-05-10 JN.1        307     0.560    500 0.614 TRUE     
#> # ℹ 35 more rows
```
