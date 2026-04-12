# Summary

`lineagefreq` is an R package for analysing pathogen lineage frequency
dynamics from genomic surveillance count data. It provides a unified
interface for fitting multinomial logistic regression models to sequence
counts, estimating growth advantages across multiple parameterisations,
generating probabilistic short-term forecasts with honest uncertainty
quantification, and evaluating forecast accuracy through rolling-origin
backtesting. The package ships with real CDC surveillance data and
supports five estimation engines — three frequentist and two Bayesian —
behind a single
[`fit_model()`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md)
interface.

# Statement of Need

The COVID-19 pandemic demonstrated that tracking the frequency
trajectories of pathogen lineages is central to public health
decision-making. When a new variant emerges, decision-makers need rapid
answers to three questions: how fast is it growing relative to
established lineages, when will it become dominant, and are current
sequencing volumes sufficient to detect it at relevant prevalence
thresholds?

These questions are answered by fitting frequency models to genomic
surveillance counts, a task that several research groups have
implemented in custom pipelines. The Nextstrain/Hubbard group developed
`evofr` for Bayesian renewal models \[@abousamra2024fitness\]. CDC
maintains internal multinomial logistic regression pipelines. Piantham,
Linton, and Nishiura described a method for converting growth rate
differences to relative reproduction numbers
\[@piantham2022predicting\]. However, no single CRAN-distributed package
has offered a unified, tested, and documented implementation covering
model fitting, multiple parameterisations, forecasting with uncertainty,
and rigorous out-of-sample evaluation.

`lineagefreq` addresses this gap. It is designed for surveillance teams,
academic epidemiology laboratories, and public health analysts who need
a self-contained, reproducible workflow with minimal setup overhead. The
package is pathogen-agnostic and works with any variant-resolved
surveillance count data.

# Key Features

## Model Fitting

The
[`fit_model()`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md)
function provides a unified interface to five estimation engines:

| Engine     | Type        | Key property                                                                             |
|------------|-------------|------------------------------------------------------------------------------------------|
| `mlr`      | Frequentist | Multinomial logistic regression via maximum likelihood                                   |
| `hier_mlr` | Frequentist | Hierarchical pooling across locations via DerSimonian–Laird random-effects meta-analysis |
| `piantham` | Frequentist | Growth rates converted to relative reproduction numbers                                  |
| `fga`      | Bayesian    | Full posterior via CmdStan                                                               |
| `garw`     | Bayesian    | Time-varying fitness via Gaussian autoregressive random walk                             |

All engines accept the same data format produced by
[`lfq_data()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md),
which validates raw count data, computes per-timepoint frequencies, and
flags unreliable periods with low sequencing volume.

## Growth Advantage Estimation

[`growth_advantage()`](https://CuiweiG.github.io/lineagefreq/reference/growth_advantage.md)
extracts fitness differences in four parameterisations: raw growth
rates, relative effective reproduction numbers (following
@piantham2022predicting), selection coefficients, and frequency doubling
times. Confidence intervals are derived from the Fisher information
matrix for frequentist engines and from posterior draws for Bayesian
engines.

## Forecasting with Uncertainty

[`forecast()`](https://CuiweiG.github.io/lineagefreq/reference/forecast.md)
generates probabilistic predictions by propagating two sources of
uncertainty: parameter estimation error (via multivariate normal
sampling from the Fisher information matrix) and multinomial sampling
variability at a user-specified effective weekly sequencing volume. The
output is a set of prediction intervals, not point forecasts, reflecting
the fundamental uncertainty in frequency trajectory projection.

## Backtesting

[`backtest()`](https://CuiweiG.github.io/lineagefreq/reference/backtest.md)
implements rolling-origin out-of-sample evaluation following the
framework described by @abousamra2024fitness. The model is repeatedly
fit on expanding historical windows, forecasts are generated at multiple
horizons, and predictions are compared to held-out observations.
[`score_forecasts()`](https://CuiweiG.github.io/lineagefreq/reference/score_forecasts.md)
computes standardised accuracy metrics (MAE, RMSE, coverage, weighted
interval score following @bracher2021evaluating), and
[`compare_models()`](https://CuiweiG.github.io/lineagefreq/reference/compare_models.md)
ranks engines by predictive performance.

## Surveillance Utilities

[`sequencing_power()`](https://CuiweiG.github.io/lineagefreq/reference/sequencing_power.md)
computes the sample size required to detect a variant at a given
prevalence with a specified confidence level.
[`summarize_emerging()`](https://CuiweiG.github.io/lineagefreq/reference/summarize_emerging.md)
scans all lineages for statistically significant growth trends using
binomial GLMs with trend tests.
[`collapse_lineages()`](https://CuiweiG.github.io/lineagefreq/reference/collapse_lineages.md)
aggregates rare lineages below a frequency or count threshold.

# Example Usage

``` r
library(lineagefreq)

# Load real CDC surveillance data
data(cdc_sarscov2_jn1)

# Prepare data
x <- lfq_data(cdc_sarscov2_jn1, lineage = lineage, date = date, count = count)

# Fit multinomial logistic regression
fit <- fit_model(x, engine = "mlr")

# Estimate relative reproduction numbers
growth_advantage(fit, type = "relative_Rt", generation_time = 5)

# Generate 4-week probabilistic forecast
fc <- forecast(fit, horizon = 28, n_sim = 1000)
autoplot(fc)

# Rolling-origin backtesting
bt <- backtest(x, engines = c("mlr", "piantham"),
               horizons = c(7, 14, 21, 28), min_train = 42)
sc <- score_forecasts(bt)
compare_models(sc)
```

# Validation

The package has been validated against published epidemiological
estimates. On CDC data covering the BA.1 to BA.2 transition, the model
recovers a relative reproduction number of approximately 1.34,
consistent with household transmission estimates from Denmark
\[@lyngse2022household\]. On the JN.1 emergence dataset, the model
reproduces the observed trajectory from below 1% to over 80% frequency.

Forecast accuracy, assessed via rolling-origin backtesting on the BA.2
transition data, yields approximately 4 percentage points MAE at a
2-week horizon and 8 percentage points at 4 weeks, with 95% prediction
interval coverage maintained across horizons.

# Datasets

The package includes four datasets: `sarscov2_us_2022` (simulated US
variant dynamics), `influenza_h3n2` (simulated H3N2 clades),
`cdc_sarscov2_jn1` (real CDC biweekly data covering the JN.1 emergence,
October 2023 to June 2024), and `cdc_ba2_transition` (real CDC data
covering the BA.1 to BA.2 replacement, December 2021 to June 2022). The
real datasets enable immediate, reproducible benchmarking against
published estimates.

# Design Philosophy

`lineagefreq` follows tidyverse conventions throughout. Broom methods
(`tidy()`, `glance()`, `augment()`) ensure integration with standard
analytical workflows. The engine registration system
([`register_engine()`](https://CuiweiG.github.io/lineagefreq/reference/register_engine.md))
allows third-party models to be plugged into the
[`fit_model()`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md)
interface, following the design philosophy of `parsnip`
\[@kuhn2022tidymodels\]. All plotting uses `ggplot2` with
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods for fitted models, forecasts, and backtest results.

# Availability

`lineagefreq` is available on CRAN at
<https://CRAN.R-project.org/package=lineagefreq> and on GitHub at
<https://github.com/CuiweiG/lineagefreq>. Four vignettes cover the core
workflow, real-data case studies, model comparison, and end-to-end
surveillance pipelines.

# References
