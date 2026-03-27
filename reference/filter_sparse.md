# Filter sparse time points and lineages

Removes time points with very low total counts and lineages observed at
too few time points.

## Usage

``` r
filter_sparse(x, min_total = 10L, min_timepoints = 3L)
```

## Arguments

- x:

  An
  [lfq_data](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
  object.

- min_total:

  Minimum total sequences per time point. Default 10.

- min_timepoints:

  Minimum number of time points a lineage must appear to be retained.
  Default 3.

## Value

An
[lfq_data](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
object with sparse entries removed.

## Examples

``` r
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
filtered <- filter_sparse(sim, min_total = 100)
```
