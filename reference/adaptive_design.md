# Adaptive sequencing allocation via Thompson sampling

Reallocates sequencing resources across regions in real time, directing
effort toward strata where uncertainty reduction has the highest
decision value. Unlike static Neyman allocation (which optimises for a
single time point), this function adapts to evolving variant dynamics
across multiple surveillance rounds.

## Usage

``` r
adaptive_design(
  data,
  capacity,
  n_rounds = NULL,
  strategy = c("thompson", "ucb"),
  target_lineage = NULL,
  exploration = 2,
  seed = NULL
)
```

## Arguments

- data:

  An `lfq_data` object with a location column.

- capacity:

  Total sequencing capacity per round (integer).

- n_rounds:

  Number of allocation rounds to simulate. Default `NULL` uses the
  number of time points in the data.

- strategy:

  Allocation strategy: `"thompson"` (default) for Thompson sampling, or
  `"ucb"` for Upper Confidence Bound.

- target_lineage:

  Character; lineage to optimise detection for. Default `NULL` optimises
  for overall frequency estimation.

- exploration:

  Exploration parameter for UCB. Default 2.0. Larger values explore
  more; smaller values exploit current best regions.

- seed:

  Random seed. Default `NULL`.

## Value

An `adaptive_allocation` S3 class with components:

- allocations:

  Tibble with `round`, `region`, `n_allocated`, `uncertainty`,
  `frequency`.

- summary:

  Tibble with per-region totals and mean allocation.

- strategy:

  Character; strategy used.

- capacity:

  Integer; per-round capacity.

## Details

**Thompson sampling** draws from the posterior distribution of frequency
estimates and allocates proportional to the posterior variance. Regions
with high uncertainty receive more sequences, but the stochastic draws
naturally balance exploration (sampling uncertain regions) and
exploitation (sampling where variants are most prevalent).

**UCB** allocates proportional to the upper confidence bound of the
estimation error: \\\text{score}\_r = \hat{\sigma}\_r + c \sqrt{2
\log(t) / n_r}\\ where \\\hat{\sigma}\_r\\ is the current estimation
uncertainty in region \\r\\, \\n_r\\ is the cumulative allocation, \\t\\
is the round, and \\c\\ is the exploration parameter.

## See also

[`surveillance_value`](https://CuiweiG.github.io/lineagefreq/reference/surveillance_value.md)
for EVOI analysis.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.3, "B" = 0.9),
  n_timepoints = 12, seed = 1)
ad <- adaptive_design(sim, capacity = 200, n_rounds = 8)
ad
#> 
#> ── Adaptive sequencing allocation 
#> Strategy: thompson | Capacity: 200/round
#> 8 rounds, 5 regions
#> 
#> # A tibble: 5 × 4
#>   region   total_allocated mean_allocated mean_uncertainty
#>   <chr>              <dbl>          <dbl>            <dbl>
#> 1 region_5             344           43            0.00105
#> 2 region_4             324           40.5          0.00105
#> 3 region_2             317           39.6          0.00105
#> 4 region_1             308           38.5          0.00105
#> 5 region_3             307           38.4          0.00105
# }
```
