# Population-level selective pressure from variant dynamics

Computes a population-level selective pressure metric from genomic
surveillance data alone, without requiring case counts or
epidemiological data. The metric quantifies how rapidly the variant
landscape is shifting and serves as an early warning signal for epidemic
growth that is robust to case underreporting.

## Usage

``` r
selective_pressure(fit, method = c("mean", "variance"))
```

## Arguments

- fit:

  An `lfq_fit` object from
  [`fit_model`](https://cuiweig.github.io/lineagefreq/reference/fit_model.md).

- method:

  Aggregation method: `"mean"` (default) for the frequency-weighted mean
  growth rate, or `"variance"` for the frequency-weighted variance of
  growth rates (higher variance indicates stronger selection).

## Value

A tibble with columns:

- date:

  Date.

- pressure:

  Selective pressure value.

- dominant_lineage:

  Lineage with highest frequency.

- dominant_freq:

  Frequency of dominant lineage.

## Details

The approach follows the framework of Figgins and Bedford (2025), where
selective pressure is defined as the variance-weighted mean growth rate
across circulating lineages: \$\$S(t) = \sum_v p_v(t) \cdot \delta_v\$\$
This represents the expected rate at which the population-average
fitness is increasing, measured entirely from sequence data.

When `method = "mean"`, the metric is positive when fitter- than-average
lineages are increasing in frequency, indicating population-level
adaptation. A sustained positive value precedes epidemic growth because
it means the effective reproduction number of the average circulating
virus is increasing.

When `method = "variance"`, the metric captures the heterogeneity of
fitness across co-circulating lineages. High variance indicates strong
directional selection; low variance indicates near-neutral drift.

This metric requires only genomic surveillance data. It does not require
case counts, hospitalisations, or wastewater data, making it applicable
in settings where epidemiological reporting is incomplete or delayed.

## References

Figgins MD, Bedford T (2025). Jointly modeling variant frequencies and
case counts to estimate relative variant severity. *medRxiv*.
[doi:10.1101/2024.12.02.24318334](https://doi.org/10.1101/2024.12.02.24318334)

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 4,
  advantages = c("A" = 1.4, "B" = 1.1, "C" = 0.8),
  n_timepoints = 15, seed = 1)
fit <- fit_model(sim, engine = "mlr")
sp <- selective_pressure(fit)
sp
#> # A tibble: 15 × 4
#>    date       pressure dominant_lineage dominant_freq
#>    <date>        <dbl> <chr>                    <dbl>
#>  1 2026-04-19    0.276 A                        0.260
#>  2 2026-04-26    0.317 A                        0.338
#>  3 2026-05-03    0.356 A                        0.421
#>  4 2026-05-10    0.392 A                        0.505
#>  5 2026-05-17    0.424 A                        0.586
#>  6 2026-05-24    0.451 A                        0.660
#>  7 2026-05-31    0.474 A                        0.725
#>  8 2026-06-07    0.492 A                        0.780
#>  9 2026-06-14    0.507 A                        0.826
#> 10 2026-06-21    0.518 A                        0.863
#> 11 2026-06-28    0.527 A                        0.893
#> 12 2026-07-05    0.534 A                        0.916
#> 13 2026-07-12    0.539 A                        0.935
#> 14 2026-07-19    0.543 A                        0.949
#> 15 2026-07-26    0.546 A                        0.961
# }
```
