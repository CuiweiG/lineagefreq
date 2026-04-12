# Validation Report

## Overview

This vignette documents the quantitative validation of `lineagefreq`
against real CDC SARS-CoV-2 surveillance data and published benchmark
results from Abousamra, Figgins, and Bedford (2024). All results are
fully reproducible using the datasets shipped with the package.

## 1. Forecast accuracy benchmark

We evaluate the MLR engine on the CDC BA.1-to-BA.2 transition dataset
(December 2021 to June 2022, 150 observations across 5 lineages) using
rolling-origin backtesting with a 42-day minimum training window.

``` r
library(lineagefreq)
data(cdc_ba2_transition)
x <- lfq_data(cdc_ba2_transition, lineage = lineage,
              date = date, count = count)
bt <- backtest(x, engines = "mlr",
               horizons = c(7, 14, 21, 28), min_train = 42)
sc <- score_forecasts(bt, metrics = c("mae", "crps", "coverage"))
```

### Results

| Horizon | Median AE | Mean AE | Coverage (95 pct) | CRPS  |
|---------|-----------|---------|-------------------|-------|
| 7 days  | 0.005     | 0.034   | 0.71              | 0.026 |
| 14 days | 0.009     | 0.057   | 0.64              | 0.040 |
| 21 days | 0.009     | 0.053   | 0.50              | 0.037 |
| 28 days | 0.020     | 0.091   | 0.40              | 0.061 |

### Comparison with published benchmarks

Abousamra et al. (2024, PLOS Computational Biology) reported MLR
forecast accuracy on a multi-country evaluation:

| Metric              | lineagefreq | Bedford Lab (2024)                          |
|---------------------|-------------|---------------------------------------------|
| Median AE at 7 days | 0.5 pp      | 0.6 pp (countries with robust surveillance) |
| Mean AE at 14 days  | 5.7 pp      | ~6 pp (approximate, from their Figure 3)    |
| Mean AE at 28 days  | 9.1 pp      | ~8-10 pp (varies by country)                |

Our point estimates are consistent with their published results. The
slightly higher MAE at 28 days may reflect the specific BA.2 transition
dynamics (four sequential subvariant sweeps occurring simultaneously),
which is a particularly challenging forecasting scenario.

## 2. Calibration diagnostics

The forecast accuracy numbers above obscure a critical issue: **the 95
percent prediction intervals are severely miscalibrated**.

``` r
cal <- calibrate(bt)
cal
plot(cal, type = "reliability")
plot(cal, type = "pit")
```

### Key finding

The PIT histogram shows a strong U-shape (KS test D = 0.41, p \< 0.001),
indicating systematic underdispersion: the MLR prediction intervals are
too narrow. At nominal 90 percent coverage, only 46 percent of
observations fall within the intervals. The mean absolute calibration
error is 0.24 — substantially worse than the 0.05 threshold that would
indicate adequate calibration.

**This finding is invisible to standard MAE-based evaluation.** A
forecaster reporting only MAE = 3.4 pp at 7 days would appear highly
accurate, but the prediction intervals — which drive decision-making
under uncertainty — are unreliable. This is precisely why calibration
diagnostics are essential.

### Recalibration

Isotonic regression recalibration widens the intervals to match
empirical coverage:

``` r
fit <- fit_model(x, engine = "mlr")
fc <- forecast(fit, horizon = 28)
fc_recal <- recalibrate(fc, bt, method = "isotonic")
```

### Conformal prediction

Conformal prediction intervals provide distribution-free coverage
guarantees. On the BA.2 data, conformal intervals are substantially
wider than parametric intervals (mean width 0.59 vs 0.07), correctly
reflecting the actual estimation uncertainty that the parametric
intervals understate.

``` r
fc_conf <- conformal_forecast(fit, x, horizon = 28,
  ci_level = 0.95, seed = 42)
```

## 3. Sequential detection

We apply SPRT-based alerting to the JN.1 emergence dataset (October 2023
to June 2024). JN.1 crossed 5 percent national frequency around November
25, 2023.

``` r
data(cdc_sarscov2_jn1)
x_jn1 <- lfq_data(cdc_sarscov2_jn1, lineage = lineage,
                   date = date, count = count)
alerts <- alert_threshold(x_jn1, method = "sprt",
  delta_1 = 0.02, alpha = 0.05)
```

### Limitation

The SPRT did not trigger an alert for JN.1 on biweekly data. This
reflects a real limitation: with only 2–3 observations between JN.1’s
first detection and its crossing of the 5 percent threshold, the
sequential test cannot accumulate sufficient evidence. Weekly or more
frequent surveillance data would provide more observations per unit time
and improve detection sensitivity. The SPRT is best suited to
high-frequency monitoring; for biweekly or monthly data, the
non-sequential
[`summarize_emerging()`](https://CuiweiG.github.io/lineagefreq/reference/summarize_emerging.md)
function using binomial trend tests is more appropriate.

## 4. Fitness decomposition

We validate the decomposition on simulated data with known ground truth.
When intrinsic transmissibility and immune escape are specified
separately, the decomposition recovers components that sum exactly to
the observed growth rate (residual \< 1e-15).

``` r
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.3, "B" = 0.9),
  n_timepoints = 20, total_per_tp = 1000, seed = 42)
fit <- fit_model(sim, engine = "mlr")
# ... construct immunity landscape and decompose ...
```

The decomposition is mathematically exact (intrinsic + escape = observed
advantage). However, the decomposition requires external immunity data —
frequency dynamics alone cannot identify the two components. This
identifiability constraint is enforced in the software:
[`fitness_decomposition()`](https://CuiweiG.github.io/lineagefreq/reference/fitness_decomposition.md)
requires an `immune_landscape` object as input and will error without
one.

## 5. Multi-pathogen applicability

The package ships with both SARS-CoV-2 (CDC surveillance) and influenza
H3N2 (simulated Nextstrain dynamics) datasets. The same `lfq_data` to
`fit_model` to `forecast` pipeline operates identically on both,
confirming pathogen-agnostic design:

``` r
data(influenza_h3n2)
x_flu <- lfq_data(influenza_h3n2, lineage = clade,
                  date = date, count = count)
fit_flu <- fit_model(x_flu, engine = "mlr")
growth_advantage(fit_flu)
```

## References

- Abousamra E, Figgins M, Bedford T (2024). Fitness models provide
  accurate short-term forecasts of SARS-CoV-2 variant frequency. *PLOS
  Computational Biology*, 20(9):e1012443.
- Lyngse FP et al. (2022). Household transmission of SARS-CoV-2 Omicron
  variant subvariants BA.1 and BA.2 in Denmark. *Nature Communications*,
  13:5760.
