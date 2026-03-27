# Plot lineage frequency model results

Plot lineage frequency model results

## Usage

``` r
# S3 method for class 'lfq_fit'
autoplot(
  object,
  type = c("frequency", "advantage", "trajectory", "residuals"),
  generation_time = NULL,
  ...
)
```

## Arguments

- object:

  An `lfq_fit` object.

- type:

  Plot type: `"frequency"` (default), `"advantage"`, `"trajectory"`, or
  `"residuals"`.

- generation_time:

  Required when `type = "advantage"`.

- ...:

  Ignored.

## Value

A ggplot object.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
fit <- fit_model(sim)
autoplot(fit)

autoplot(fit, type = "advantage", generation_time = 5)

# }
```
