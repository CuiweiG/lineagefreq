# Recalibrate prediction intervals

Applies post-hoc recalibration to improve the coverage properties of
prediction intervals from
[`forecast`](https://CuiweiG.github.io/lineagefreq/reference/forecast.md).
Two methods are available: isotonic regression (nonparametric,
monotonicity- preserving) and Platt scaling (logistic, parametric).

## Usage

``` r
recalibrate(forecast_obj, bt, method = c("isotonic", "platt"))
```

## Arguments

- forecast_obj:

  An `lfq_forecast` object to recalibrate.

- bt:

  An `lfq_backtest` object providing the calibration data. The mapping
  from nominal to empirical coverage is learned from the backtest
  residuals.

- method:

  Recalibration method: `"isotonic"` (default) for isotonic regression
  on the empirical coverage function, or `"platt"` for Platt scaling via
  logistic regression.

## Value

An `lfq_forecast` object with recalibrated `.lower` and `.upper` bounds.
The object retains all attributes of the original forecast and can be
passed to
[`autoplot`](https://ggplot2.tidyverse.org/reference/autoplot.html).

## Details

**Isotonic regression:** Learns the monotone mapping from nominal
coverage levels to observed coverage using the backtest data, then
inverts it to find the nominal level that achieves the desired empirical
coverage. This is a nonparametric approach that requires no
distributional assumptions.

**Platt scaling:** Fits a logistic regression of the form
\\P(\text{covered}) = \text{logit}^{-1}(a \cdot z + b)\\ where \\z\\ is
the standardised residual, and uses the fitted model to adjust
prediction interval widths.

Both methods require a backtest object with sufficient data to estimate
the calibration mapping reliably (at least 30 forecast-observation pairs
recommended).

## References

Platt JC (1999). Probabilistic outputs for support vector machines and
comparisons to regularized likelihood methods. *Advances in Large Margin
Classifiers*, 61–74.

## See also

[`calibrate`](https://CuiweiG.github.io/lineagefreq/reference/calibrate.md)
for calibration diagnostics.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.2, "B" = 0.8),
  n_timepoints = 20, seed = 1)
fit <- fit_model(sim, engine = "mlr")
fc  <- forecast(fit, horizon = 21)
bt  <- backtest(sim, engines = "mlr",
  horizons = c(7, 14), min_train = 42)
fc_recal <- recalibrate(fc, bt, method = "isotonic")
# }
```
