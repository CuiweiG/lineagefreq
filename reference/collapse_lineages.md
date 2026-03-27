# Collapse rare lineages into an aggregate group

Merges lineages that never exceed a frequency or count threshold into a
single group (default "Other"). Useful for reducing noise from dozens of
low-frequency lineages.

## Usage

``` r
collapse_lineages(
  x,
  min_freq = 0.01,
  min_count = 10L,
  other_label = "Other",
  mapping = NULL
)
```

## Arguments

- x:

  An
  [lfq_data](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
  object.

- min_freq:

  Minimum peak frequency a lineage must reach at any time point to be
  kept. Default 0.01 (1%).

- min_count:

  Minimum total count across all time points to be kept. Default 10.

- other_label:

  Label for the collapsed group. Default "Other".

- mapping:

  Optional named character vector for custom grouping. Names = original
  lineage names, values = group names. If provided, `min_freq` and
  `min_count` are ignored.

## Value

An
[lfq_data](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
object with rare lineages merged.

## Examples

``` r
sim <- simulate_dynamics(n_lineages = 6,
  advantages = c(A = 1.3, B = 1.1, C = 0.95, D = 0.5, E = 0.3),
  n_timepoints = 12, seed = 1)
collapsed <- collapse_lineages(sim, min_freq = 0.05)
attr(collapsed, "lineages")
#> [1] "A"   "B"   "C"   "D"   "E"   "ref"
```
