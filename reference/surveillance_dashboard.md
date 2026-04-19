# Comprehensive surveillance quality dashboard

Produces a multi-panel display combining calibration diagnostics,
detection power, estimation quality, and current variant landscape into
a single figure suitable for weekly surveillance reports. Designed for
programme managers rather than statisticians.

## Usage

``` r
surveillance_dashboard(fit, data, bt = NULL, target_prevalence = 0.01)
```

## Arguments

- fit:

  An `lfq_fit` object from
  [`fit_model`](https://cuiweig.github.io/lineagefreq/reference/fit_model.md).

- data:

  The `lfq_data` object used for fitting.

- bt:

  Optional `lfq_backtest` object for calibration panel. If `NULL`, the
  calibration panel is omitted.

- target_prevalence:

  Prevalence for detection power calculation. Default 0.01 (1 percent).

## Value

A list of ggplot objects with class `surveillance_dashboard`. A print
method renders all panels.

## Details

The dashboard contains up to four panels: (1) current frequency
landscape, (2) growth advantage forest plot, (3) detection power curve,
and (4) calibration reliability diagram (if backtest data are provided).

## See also

[`surveillance_value`](https://cuiweig.github.io/lineagefreq/reference/surveillance_value.md)
for EVOI analysis,
[`alert_threshold`](https://cuiweig.github.io/lineagefreq/reference/alert_threshold.md)
for sequential detection.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.3, "B" = 0.9),
  n_timepoints = 15, seed = 1)
fit <- fit_model(sim, engine = "mlr")
panels <- surveillance_dashboard(fit, sim)
# }
```
