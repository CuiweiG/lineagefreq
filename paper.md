# Summary

`lineagefreq` is an R package for analysing the frequency dynamics of
pathogen lineages from genomic surveillance count data. It provides a
unified interface for fitting multinomial logistic regression models to
sequence counts, estimating lineage-specific growth advantages in
multiple parameterisations, generating probabilistic short-term
frequency forecasts with honest uncertainty quantification, and
evaluating forecast accuracy through rolling-origin backtesting. The
package supports five estimation engines — three frequentist and two
Bayesian — behind a single
[`fit_model()`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md)
function, and ships with real U.S. CDC SARS-CoV-2 surveillance data for
immediate reproducibility.

# Statement of Need

When a novel pathogen variant emerges, public health authorities require
rapid answers to three operational questions: how fast is the new
lineage growing relative to established variants, when will it become
the dominant circulating strain, and are current sequencing volumes
adequate to detect it at decision-relevant prevalence thresholds? These
questions arise with every variant of concern — from SARS-CoV-2 Omicron
sublineages to seasonal influenza clade replacements — and they are
fundamentally questions about lineage frequency dynamics that can be
addressed with well-characterised statistical models.

The multinomial logistic regression (MLR) framework provides a natural
model for this setting. Observed sequence counts at each surveillance
time point are treated as draws from a multinomial distribution, with
log-linear frequency trajectories parameterised by lineage-specific
intercepts and growth rates. The growth rate difference between a focal
lineage and a reference lineage can be converted to a relative effective
reproduction number using the generation time
\[@piantham2022predicting\]. This approach has been validated
extensively during the COVID-19 pandemic, with several groups
implementing custom pipelines for real-time variant tracking
\[@abousamra2024fitness\].

Despite the maturity of the underlying methodology, no CRAN-distributed
R package has previously offered a unified, documented implementation
covering model fitting across multiple engines, growth advantage
estimation in multiple parameterisations, probabilistic forecasting with
rigorous uncertainty propagation, and standardised out-of-sample
evaluation. Existing tools are either embedded in specialised Bayesian
frameworks requiring Stan infrastructure \[@carpenter2017stan\],
maintained as internal pipelines within public health agencies, or
implemented as single-purpose scripts that are difficult to reproduce
across institutions.

`lineagefreq` fills this gap. It is designed for surveillance teams,
academic epidemiology laboratories, and public health analysts who need
a self-contained, tested, and documented solution with minimal setup
requirements. The package is pathogen-agnostic and handles any
variant-resolved surveillance count data.

# Methodology

## Model Specification

The core model treats sequence counts $Y_{v}(t)$ for lineage $v$ at time
$t$ as multinomially distributed with frequency parameters:

$$p_{v}(t) = \frac{\exp\left( \alpha_{v} + \delta_{v} \cdot t \right)}{\sum\limits_{k}\exp\left( \alpha_{k} + \delta_{k} \cdot t \right)}$$

where $\alpha_{v}$ is an intercept governing initial frequency and
$\delta_{v}$ is the growth rate per unit time. A designated pivot
lineage has $\alpha_{\text{pivot}} = \delta_{\text{pivot}} = 0$, serving
as the reference against which all growth advantages are measured.
Parameters are estimated by maximum likelihood via BFGS optimisation,
with confidence intervals derived from the Fisher information matrix.

## Estimation Engines

Five engines are available behind the
[`fit_model()`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md)
interface:

- **`mlr`**: Standard multinomial logistic regression via maximum
  likelihood. Suitable for single-location analyses with adequate sample
  sizes.
- **`hier_mlr`**: Hierarchical pooling across multiple locations using
  DerSimonian–Laird random-effects meta-analysis
  \[@dersimonian1986meta\]. Fits MLR independently per location, then
  shrinks location-specific growth rate estimates toward the pooled
  mean. Particularly valuable during early-phase surveillance when
  individual jurisdictions have few sequences.
- **`piantham`**: Wraps the MLR engine and converts growth rates to
  relative effective reproduction numbers following
  @piantham2022predicting, requiring specification of the generation
  time.
- **`fga`** and **`garw`**: Bayesian engines using Stan
  \[@carpenter2017stan\] via `cmdstanr`, providing full posterior
  distributions. The `garw` (Gaussian autoregressive random walk) engine
  additionally accommodates time-varying fitness, relaxing the
  constant-growth-rate assumption.

The engine registration system
([`register_engine()`](https://CuiweiG.github.io/lineagefreq/reference/register_engine.md),
[`unregister_engine()`](https://CuiweiG.github.io/lineagefreq/reference/unregister_engine.md))
permits third-party models to be integrated into the same
[`fit_model()`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md)
interface, following the extensible design philosophy of `parsnip`
\[@kuhn2022tidymodels\].

## Growth Advantage Estimation

The
[`growth_advantage()`](https://CuiweiG.github.io/lineagefreq/reference/growth_advantage.md)
function extracts fitness differences in four parameterisations:

- **Growth rate**: the raw $\delta_{v}$ per time-scale unit.
- **Relative reproduction number**:
  $R_{t}^{(v)} = \exp\left( \delta_{v} \cdot \tau/T \right)$, where
  $\tau$ is the generation time and $T$ is the time scale
  \[@piantham2022predicting\].
- **Selection coefficient**: $R_{t}^{(v)} - 1$, expressing the fitness
  advantage as a proportion.
- **Doubling time**: the number of days for the frequency ratio relative
  to the pivot to double.

Confidence intervals use the delta method applied to the
variance-covariance matrix from maximum likelihood estimation, computed
in log-space for ratio quantities to improve accuracy.

## Forecasting

The
[`forecast()`](https://CuiweiG.github.io/lineagefreq/reference/forecast.md)
function generates probabilistic predictions by propagating two sources
of uncertainty. Parameter estimation error is captured by drawing from a
multivariate normal distribution centred on the maximum likelihood
estimates with covariance given by the inverse Fisher information.
Multinomial sampling variability is optionally added at a user-specified
effective weekly sequencing volume. The result is a set of prediction
intervals computed from quantiles across simulation draws — not point
forecasts.

## Backtesting

The
[`backtest()`](https://CuiweiG.github.io/lineagefreq/reference/backtest.md)
function implements rolling-origin out-of-sample evaluation following
the framework described by @abousamra2024fitness. At each origin date,
the model is fit on historical data and forecasts are generated at
multiple horizons (default: 7, 14, 21, and 28 days). Predictions are
compared to held-out observations using
[`score_forecasts()`](https://CuiweiG.github.io/lineagefreq/reference/score_forecasts.md),
which computes standardised accuracy metrics: mean absolute error, root
mean squared error, prediction interval coverage, and weighted interval
score \[@bracher2021evaluating\]. The
[`compare_models()`](https://CuiweiG.github.io/lineagefreq/reference/compare_models.md)
function ranks engines by predictive performance.

# Surveillance Utilities

Two additional functions address practical surveillance operations.
[`sequencing_power()`](https://CuiweiG.github.io/lineagefreq/reference/sequencing_power.md)
computes the sample size required to detect a lineage at a specified
prevalence with a given confidence level, answering the recurring
question of whether a jurisdiction’s sequencing volume is adequate.
[`summarize_emerging()`](https://CuiweiG.github.io/lineagefreq/reference/summarize_emerging.md)
scans all lineages for statistically significant growth trends using
binomial GLMs with time-trend coefficients, flagging those exceeding a
user-specified frequency threshold. Preprocessing functions
[`collapse_lineages()`](https://CuiweiG.github.io/lineagefreq/reference/collapse_lineages.md)
and
[`filter_sparse()`](https://CuiweiG.github.io/lineagefreq/reference/filter_sparse.md)
handle low-frequency lineage aggregation and removal of time points with
insufficient sequencing depth.

# Included Datasets

The package ships with four datasets enabling immediate, reproducible
analysis:

- **`cdc_sarscov2_jn1`**: Real CDC biweekly surveillance data covering
  the JN.1 emergence (October 2023 to June 2024), with 10 lineages
  including JN.1, XBB.1.5, and EG.5.1. Source: CDC COVID Data Tracker,
  public domain.
- **`cdc_ba2_transition`**: Real CDC data covering the
  well-characterised BA.1 to BA.2 replacement (December 2021 to June
  2022), providing an independent validation dataset with published
  growth advantage estimates for comparison.
- **`sarscov2_us_2022`**: Simulated weekly U.S. data representing four
  sequential variant sweeps (BA.1, BA.2, BA.4/5, BQ.1).
- **`influenza_h3n2`**: Simulated seasonal influenza H3N2 clade data
  based on Nextstrain dynamics.

# Validation

The package has been validated against published epidemiological
estimates. On CDC data covering the BA.1 to BA.2 transition, the model
recovers a relative reproduction number of approximately 1.34,
consistent with household transmission estimates of 1.3–1.5 from Denmark
\[@lyngse2022household\]. On the JN.1 emergence dataset, the model
reproduces the observed frequency trajectory from below 1% to over 80%.

Forecast accuracy, assessed through rolling-origin backtesting on the
BA.2 transition data, yields approximately 4 percentage points MAE at a
2-week horizon and 8 percentage points at 4 weeks. The 95% prediction
interval coverage is maintained across forecast horizons, indicating
well-calibrated uncertainty quantification.

# Example Usage

``` r
library(lineagefreq)

# Load real CDC surveillance data
data(cdc_sarscov2_jn1)

# Prepare and validate count data
x <- lfq_data(cdc_sarscov2_jn1, lineage = lineage,
              date = date, count = count)

# Fit multinomial logistic regression
fit <- fit_model(x, engine = "mlr")

# Estimate relative reproduction numbers
growth_advantage(fit, type = "relative_Rt", generation_time = 5)

# Generate 4-week probabilistic forecast
fc <- forecast(fit, horizon = 28, n_sim = 1000)
autoplot(fc)

# Rolling-origin backtesting with standardised scoring
bt <- backtest(x, engines = c("mlr", "piantham"),
               horizons = c(7, 14, 21, 28), min_train = 42)
sc <- score_forecasts(bt)
compare_models(sc)
```

# Software Design

`lineagefreq` follows tidyverse conventions throughout. Broom methods
(`tidy()`, `glance()`, `augment()`) ensure that model outputs integrate
into standard analytical pipelines. All plotting uses `ggplot2` with
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods for fitted models, forecasts, and backtest summaries. The
package uses tibble subclassing for data containers and rlang for tidy
evaluation.

The package is available on CRAN at
<https://CRAN.R-project.org/package=lineagefreq> and on GitHub at
<https://github.com/CuiweiG/lineagefreq>. Four vignettes cover the core
workflow, real-data case studies, model comparison, and end-to-end
surveillance pipelines. The test suite comprises 251 unit tests.

# Acknowledgements

The methodological foundation draws on work by Abousamra, Figgins, and
Bedford on fitness model evaluation, and by Piantham, Linton, and
Nishiura on relative reproduction number approximation. The CDC COVID
Data Tracker provided the real surveillance data included in the
package.

# References
