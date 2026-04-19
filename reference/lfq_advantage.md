# Pipe-friendly growth advantage extraction

Pipe-friendly growth advantage extraction

## Usage

``` r
lfq_advantage(fit, type = "relative_Rt", generation_time = NULL, ...)
```

## Arguments

- fit:

  An `lfq_fit` object.

- type:

  Output type. Default `"relative_Rt"`.

- generation_time:

  Mean generation time in days.

- ...:

  Passed to
  [`growth_advantage()`](https://cuiweig.github.io/lineagefreq/reference/growth_advantage.md).

## Value

A tibble of growth advantages.
