# Pseudo-Prospective Evaluation of Conformal Prediction

Simulates real-time deployment: at each origin, all preceding origins
form the conformal calibration set, and the forecast at the next origin
is evaluated. No future data is ever used for calibration — this is a
genuine expanding-window evaluation.

## Usage

``` r
evaluate_prospective(
  data,
  engine = "mlr",
  horizons = c(7L, 14L),
  min_train = 42L,
  min_cal = 5L,
  ci_level = 0.95,
  gamma = 0.05,
  ...
)
```

## Arguments

- data:

  An `lfq_data` object.

- engine:

  Character; estimation engine (default `"mlr"`).

- horizons:

  Integer vector; forecast horizons in days (default `c(7L, 14L)`).

- min_train:

  Integer; minimum training window in days (default 42).

- min_cal:

  Integer; minimum calibration origins before conformal intervals are
  computed (default 5).

- ci_level:

  Numeric in (0,1); nominal coverage (default 0.95).

- gamma:

  Numeric; ACI learning rate (default 0.05).

- ...:

  Passed to
  [`fit_model()`](https://cuiweig.github.io/lineagefreq/reference/fit_model.md).

## Value

An `lfq_prospective` S3 object (list) with:

- results:

  Tibble with origin_date, target_date, horizon, lineage, predicted,
  observed, lower_param, upper_param, lower_static, upper_static,
  lower_aci, upper_aci, radius_static, radius_aci, alpha_aci.

- coverage:

  Tibble summarising cumulative coverage by method (parametric, static
  conformal, ACI) at each step.

- summary:

  Tibble with method, overall coverage, mean width, Winkler score.

- n_origins:

  Total number of evaluation origins.

## Details

At each origin \\t_i\\ with \\i \geq\\ `min_cal`:

1.  Fit the model on data up to \\t_i\\.

2.  Compute residuals at all previous origins \\t_1, \ldots, t\_{i-1}\\
    (the calibration set).

3.  Static conformal: radius = \\(1-\alpha)(1+1/n)\\ quantile of
    absolute calibration residuals.

4.  ACI: radius adjusted by online update of \\\alpha_t\\.

5.  Evaluate coverage at \\t\_{i+1}\\.

This differs from
[`conformal_forecast()`](https://cuiweig.github.io/lineagefreq/reference/conformal_forecast.md)
which uses a fixed temporal split. Here, the calibration set expands
over time, mimicking a surveillance agency accumulating validation data.
