# Conformal prediction intervals for lineage frequencies

Produces distribution-free prediction intervals with finite-sample
coverage guarantees using split conformal inference. Unlike the
parametric intervals from
[`forecast`](https://CuiweiG.github.io/lineagefreq/reference/forecast.md),
conformal intervals require no distributional assumptions on the
residuals and are valid under exchangeability.

## Usage

``` r
conformal_forecast(
  fit,
  data,
  horizon = 28L,
  ci_level = 0.95,
  method = c("split", "aci"),
  cal_fraction = 0.3,
  gamma = 0.05,
  seed = NULL
)
```

## Arguments

- fit:

  An `lfq_fit` object (any engine).

- data:

  The `lfq_data` object used to fit the model.

- horizon:

  Number of days to forecast. Default 28.

- ci_level:

  Target coverage level. Default 0.95.

- method:

  Conformal method: `"split"` (default) for split conformal prediction,
  or `"aci"` for adaptive conformal inference with online coverage
  correction.

- cal_fraction:

  Fraction of the data reserved for the calibration set (split conformal
  only). Default 0.3.

- gamma:

  Learning rate for adaptive conformal inference. Default 0.05. Controls
  how quickly the coverage target adjusts in response to observed
  miscoverage.

- seed:

  Random seed for the calibration split. Default `NULL`.

## Value

An `lfq_forecast` object with conformal prediction intervals. The object
is fully compatible with
[`autoplot`](https://ggplot2.tidyverse.org/reference/autoplot.html) and
other forecast methods.

## Details

**Split conformal prediction** partitions the training data into a
proper training set and a calibration set. The model is refit on the
training set, and conformity scores (absolute residuals) are computed on
the calibration set. The prediction interval at a new point is the point
forecast plus or minus the \\(1 - \alpha)(1 + 1/n\_{\text{cal}})\\
quantile of the calibration scores. This provides exact \\1 - \alpha\\
marginal coverage under exchangeability (Vovk et al. 2005).

**Adaptive conformal inference (ACI)** (Gibbs and Candes, 2021) adjusts
the miscoverage level online to maintain long-run coverage even when the
data distribution shifts over time, as is typical in surveillance data
during variant replacement waves.

## References

Vovk V, Gammerman A, Shafer G (2005). *Algorithmic Learning in a Random
World*. Springer.

Gibbs I, Candes E (2021). Adaptive conformal inference under
distribution shift. *Advances in Neural Information Processing Systems*,
34.

## See also

[`forecast`](https://CuiweiG.github.io/lineagefreq/reference/forecast.md)
for parametric prediction intervals,
[`calibrate`](https://CuiweiG.github.io/lineagefreq/reference/calibrate.md)
for calibration diagnostics.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8),
  n_timepoints = 20, seed = 1)
fit <- fit_model(sim, engine = "mlr")
fc_conf <- conformal_forecast(fit, sim, horizon = 21)
fc_conf
#> 
#> ── Lineage frequency forecast 
#> Engine: mlr
#> Forecast start: 2026-08-24 | Horizon: 21 days
#> CI level: 95%
#> 60 fitted + 9 forecast rows
#> 
#> # A tibble: 69 × 6
#>    .date      .lineage .median .lower .upper .type 
#>    <date>     <chr>      <dbl>  <dbl>  <dbl> <chr> 
#>  1 2026-04-12 A          0.341  0.341  0.341 fitted
#>  2 2026-04-12 B          0.324  0.324  0.324 fitted
#>  3 2026-04-12 ref        0.335  0.335  0.335 fitted
#>  4 2026-04-19 A          0.408  0.408  0.408 fitted
#>  5 2026-04-19 B          0.257  0.257  0.257 fitted
#>  6 2026-04-19 ref        0.335  0.335  0.335 fitted
#>  7 2026-04-26 A          0.476  0.476  0.476 fitted
#>  8 2026-04-26 B          0.198  0.198  0.198 fitted
#>  9 2026-04-26 ref        0.326  0.326  0.326 fitted
#> 10 2026-05-03 A          0.541  0.541  0.541 fitted
#> # ℹ 59 more rows
# }
```
