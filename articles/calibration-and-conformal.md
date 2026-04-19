# Prediction Calibration and Conformal Inference

## Why calibration matters

A prediction interval that claims 95% coverage but actually covers the
true value only 70% of the time is worse than useless — it gives
decision-makers false confidence. In genomic surveillance, this can mean
underestimating the probability that a variant reaches a critical
threshold, leading to delayed public health response.

Calibration asks a simple question: when a model says “there is a 95%
chance the true frequency lies in this interval,” does the observation
actually fall within that interval 95% of the time? Surprisingly few
surveillance forecasting tools assess this systematically.

`lineagefreq` provides three capabilities for forecast calibration that,
to our knowledge, no other genomic surveillance package offers:
calibration diagnostics
([`calibrate()`](https://cuiweig.github.io/lineagefreq/reference/calibrate.md)),
post-hoc recalibration
([`recalibrate()`](https://cuiweig.github.io/lineagefreq/reference/recalibrate.md)),
and distribution-free conformal prediction intervals
([`conformal_forecast()`](https://cuiweig.github.io/lineagefreq/reference/conformal_forecast.md)).

## PIT diagnostics on real CDC data

The Probability Integral Transform (PIT) provides a comprehensive
calibration diagnostic. For a perfectly calibrated forecaster, PIT
values are uniformly distributed on \[0, 1\]. Departures from uniformity
reveal specific miscalibration patterns: a U-shaped PIT histogram
indicates underdispersion (intervals too narrow), while a hump shape
indicates overdispersion (intervals too wide).

``` r
library(lineagefreq)
data(cdc_ba2_transition)

x <- lfq_data(cdc_ba2_transition, lineage = lineage,
              date = date, count = count)
bt <- backtest(x, engines = "mlr",
               horizons = c(7, 14, 21, 28), min_train = 42)

cal <- calibrate(bt)
cal

# PIT histogram
plot(cal, type = "pit")

# Reliability diagram
plot(cal, type = "reliability")
```

The reliability diagram plots observed coverage against nominal coverage
at levels 10% through 90%. Perfect calibration lies on the diagonal.
Points above the diagonal indicate overconfidence (the model claims
narrower intervals than warranted); points below indicate conservatism.

## Recalibration

When the PIT diagnostic reveals miscalibration,
[`recalibrate()`](https://cuiweig.github.io/lineagefreq/reference/recalibrate.md)
adjusts prediction intervals to improve coverage. Two methods are
available:

**Isotonic regression** learns a nonparametric, monotone mapping from
nominal to observed coverage using the backtest data and inverts it.
This requires no distributional assumptions.

**Platt scaling** fits a logistic regression to the relationship between
standardised residuals and coverage, providing a smooth parametric
adjustment.

``` r
fit <- fit_model(x, engine = "mlr")
fc  <- forecast(fit, horizon = 28)

# Isotonic recalibration
fc_recal <- recalibrate(fc, bt, method = "isotonic")
autoplot(fc_recal)
```

## Conformal prediction intervals

Parametric prediction intervals from
[`forecast()`](https://cuiweig.github.io/lineagefreq/reference/forecast.md)
assume that the multinomial logistic model is correctly specified and
that the parameter estimation uncertainty is well captured by the Fisher
information matrix. When these assumptions are violated — as they may be
during rapid variant replacement waves — the intervals can be
substantially miscalibrated.

Conformal prediction offers a distribution-free alternative. Split
conformal inference (Vovk et al. 2005) provides prediction intervals
with exact finite-sample coverage guarantees under the assumption of
exchangeability:

``` r
fc_conf <- conformal_forecast(fit, x, horizon = 28,
                              ci_level = 0.95, seed = 42)
autoplot(fc_conf)
```

The conformal intervals are typically wider than the parametric
intervals when the model is well-specified, but narrower when the model
is misspecified — precisely the situation where accurate uncertainty
quantification matters most.

### Adaptive conformal inference

For time-series data where the distribution shifts over time (as during
a variant replacement wave), the exchangeability assumption is violated.
Adaptive conformal inference (ACI; Gibbs and Candes, 2021) addresses
this by adjusting the miscoverage level online:

``` r
fc_aci <- conformal_forecast(fit, x, horizon = 28,
                             method = "aci", gamma = 0.05)
```

The `gamma` parameter controls the learning rate: larger values track
distribution shifts more aggressively but increase interval variability.

## Proper scoring rules

[`score_forecasts()`](https://cuiweig.github.io/lineagefreq/reference/score_forecasts.md)
now supports four additional proper scoring rules alongside the original
MAE, RMSE, coverage, and WIS:

- **CRPS** (Continuous Ranked Probability Score): the integral of the
  squared difference between the forecast CDF and the empirical CDF of
  the observation. Simultaneously rewards calibration and sharpness.
- **Log score**: the negative log-likelihood of the observation under
  the forecast distribution. Heavily penalises confident wrong
  forecasts.
- **DSS** (Dawid-Sebastiani Score): depends only on the first two
  moments, providing a computationally efficient alternative to the log
  score.
- **Calibration score**: mean squared calibration error across nominal
  levels 10%–90%.

``` r
sc <- score_forecasts(bt, metrics = c("mae", "crps", "log_score",
                                      "dss", "calibration"))
compare_models(sc)
```

These scoring rules are proper in the sense of Gneiting and Raftery
(2007): they are minimised in expectation when the forecast distribution
equals the true data-generating distribution. This means that a
forecaster cannot improve their expected score by issuing dishonest
predictions.

## References

- Gneiting T, Balabdaoui F, Raftery AE (2007). Probabilistic forecasts,
  calibration and sharpness. *JRSS-B* 69(2):243–268.
- Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
  prediction, and estimation. *JASA* 102(477):359–378.
- Vovk V, Gammerman A, Shafer G (2005). *Algorithmic Learning in a
  Random World*. Springer.
- Gibbs I, Candes E (2021). Adaptive conformal inference under
  distribution shift. *NeurIPS* 34.
