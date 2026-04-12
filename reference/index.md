# Package index

## Data

Create and manipulate lineage frequency data

- [`lfq_data()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
  : Create a lineage frequency data object
- [`is_lfq_data()`](https://CuiweiG.github.io/lineagefreq/reference/is_lfq_data.md)
  : Test if an object is an lfq_data object
- [`as_lfq_data()`](https://CuiweiG.github.io/lineagefreq/reference/as_lfq_data.md)
  : Coerce to lfq_data
- [`as.data.frame(`*`<lfq_data>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/as.data.frame.lfq_data.md)
  : Convert lfq_data to long-format tibble
- [`read_lineage_counts()`](https://CuiweiG.github.io/lineagefreq/reference/read_lineage_counts.md)
  : Read lineage count data from a CSV file
- [`simulate_dynamics()`](https://CuiweiG.github.io/lineagefreq/reference/simulate_dynamics.md)
  : Simulate lineage frequency dynamics
- [`collapse_lineages()`](https://CuiweiG.github.io/lineagefreq/reference/collapse_lineages.md)
  : Collapse rare lineages into an aggregate group
- [`filter_sparse()`](https://CuiweiG.github.io/lineagefreq/reference/filter_sparse.md)
  : Filter sparse time points and lineages

## Modeling

Fit frequency dynamics models

- [`fit_model()`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md)
  : Fit a lineage frequency model
- [`lfq_engines()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_engines.md)
  : List available modeling engines
- [`register_engine()`](https://CuiweiG.github.io/lineagefreq/reference/register_engine.md)
  : Register a custom modeling engine
- [`unregister_engine()`](https://CuiweiG.github.io/lineagefreq/reference/unregister_engine.md)
  : Remove a registered engine

## Pipe API

Tidyverse-style chaining

- [`lfq_fit()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_fit.md)
  : Pipe-friendly model fitting
- [`lfq_advantage()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_advantage.md)
  : Pipe-friendly growth advantage extraction
- [`lfq_forecast()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_forecast.md)
  : Pipe-friendly forecasting
- [`lfq_score()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_score.md)
  : Pipe-friendly backtesting + scoring
- [`lfq_summary()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_summary.md)
  : Convert lfq_fit results to a summary tibble

## Inference

Extract growth advantages and forecasts

- [`growth_advantage()`](https://CuiweiG.github.io/lineagefreq/reference/growth_advantage.md)
  : Extract growth advantage estimates
- [`forecast()`](https://CuiweiG.github.io/lineagefreq/reference/forecast.md)
  : Forecast lineage frequencies (generic)
- [`forecast(`*`<lfq_fit>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/forecast.lfq_fit.md)
  : Forecast lineage frequencies
- [`summarize_emerging()`](https://CuiweiG.github.io/lineagefreq/reference/summarize_emerging.md)
  : Summarize emerging lineages
- [`sequencing_power()`](https://CuiweiG.github.io/lineagefreq/reference/sequencing_power.md)
  : Sequencing power analysis

## Calibration & Uncertainty

Assess and correct prediction interval calibration

- [`calibrate()`](https://CuiweiG.github.io/lineagefreq/reference/calibrate.md)
  : Calibration diagnostics for lineage frequency forecasts
- [`conformal_forecast()`](https://CuiweiG.github.io/lineagefreq/reference/conformal_forecast.md)
  : Conformal prediction intervals for lineage frequencies
- [`recalibrate()`](https://CuiweiG.github.io/lineagefreq/reference/recalibrate.md)
  : Recalibrate prediction intervals
- [`plot(`*`<calibration_report>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/plot.calibration_report.md)
  : Plot calibration diagnostics

## Fitness & Immunity

Decompose growth advantage into transmissibility and immune escape

- [`fitness_decomposition()`](https://CuiweiG.github.io/lineagefreq/reference/fitness_decomposition.md)
  : Decompose variant fitness into transmissibility and immune escape
- [`immune_landscape()`](https://CuiweiG.github.io/lineagefreq/reference/immune_landscape.md)
  : Construct a population immunity landscape
- [`fit_dms_prior()`](https://CuiweiG.github.io/lineagefreq/reference/fit_dms_prior.md)
  : Fit model with Deep Mutational Scanning priors
- [`selective_pressure()`](https://CuiweiG.github.io/lineagefreq/reference/selective_pressure.md)
  : Population-level selective pressure from variant dynamics
- [`plot(`*`<fitness_decomposition>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/plot.fitness_decomposition.md)
  : Plot fitness decomposition
- [`plot(`*`<immune_landscape>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/plot.immune_landscape.md)
  : Plot population immunity landscape
- [`tidy.fitness_decomposition()`](https://CuiweiG.github.io/lineagefreq/reference/tidy.fitness_decomposition.md)
  : Tidy a fitness decomposition

## Surveillance Optimization

Information-theoretic resource allocation and emergence detection

- [`surveillance_value()`](https://CuiweiG.github.io/lineagefreq/reference/surveillance_value.md)
  : Expected Value of Information for genomic surveillance
- [`adaptive_design()`](https://CuiweiG.github.io/lineagefreq/reference/adaptive_design.md)
  : Adaptive sequencing allocation via Thompson sampling
- [`alert_threshold()`](https://CuiweiG.github.io/lineagefreq/reference/alert_threshold.md)
  : Sequential detection of emerging variants
- [`detection_horizon()`](https://CuiweiG.github.io/lineagefreq/reference/detection_horizon.md)
  : Detection horizon for an emerging variant
- [`surveillance_dashboard()`](https://CuiweiG.github.io/lineagefreq/reference/surveillance_dashboard.md)
  : Comprehensive surveillance quality dashboard
- [`plot(`*`<adaptive_allocation>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/plot.adaptive_allocation.md)
  : Plot adaptive allocation
- [`plot(`*`<evoi>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/plot.evoi.md)
  : Plot EVOI curve

## Backtesting

Evaluate forecast accuracy

- [`backtest()`](https://CuiweiG.github.io/lineagefreq/reference/backtest.md)
  : Rolling-origin backtesting of lineage frequency models
- [`score_forecasts()`](https://CuiweiG.github.io/lineagefreq/reference/score_forecasts.md)
  : Score backtest forecast accuracy
- [`compare_models()`](https://CuiweiG.github.io/lineagefreq/reference/compare_models.md)
  : Compare model engines from backtest scores

## Visualization

Publication-ready plots

- [`autoplot(`*`<lfq_fit>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/autoplot.lfq_fit.md)
  : Plot lineage frequency model results
- [`autoplot(`*`<lfq_forecast>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/autoplot.lfq_forecast.md)
  : Plot a lineage frequency forecast
- [`plot_backtest()`](https://CuiweiG.github.io/lineagefreq/reference/plot_backtest.md)
  : Plot backtest scores

## S3 Methods

Standard R model interface

- [`print(`*`<lfq_fit>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/print.lfq_fit.md)
  : Print a lineage frequency model
- [`summary(`*`<lfq_fit>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/summary.lfq_fit.md)
  : Summarise a lineage frequency model
- [`coef(`*`<lfq_fit>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/coef.lfq_fit.md)
  : Extract coefficients from a lineage frequency model
- [`tidy.lfq_fit()`](https://CuiweiG.github.io/lineagefreq/reference/tidy.lfq_fit.md)
  : Tidy an lfq_fit object
- [`glance.lfq_fit()`](https://CuiweiG.github.io/lineagefreq/reference/glance.lfq_fit.md)
  : Glance at an lfq_fit object
- [`augment.lfq_fit()`](https://CuiweiG.github.io/lineagefreq/reference/augment.lfq_fit.md)
  : Augment data with fitted values from an lfq_fit object

## Utilities

- [`lfq_stan_available()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_stan_available.md)
  : Check if 'CmdStan' backend is available
- [`lfq_version()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_version.md)
  : Package version and system information

## Datasets

Built-in example datasets

- [`sarscov2_us_2022`](https://CuiweiG.github.io/lineagefreq/reference/sarscov2_us_2022.md)
  : Simulated SARS-CoV-2 variant frequency data (US, 2022)
- [`cdc_ba2_transition`](https://CuiweiG.github.io/lineagefreq/reference/cdc_ba2_transition.md)
  : CDC SARS-CoV-2 variant proportions: BA.1 to BA.2 transition (US,
  2022)
- [`cdc_sarscov2_jn1`](https://CuiweiG.github.io/lineagefreq/reference/cdc_sarscov2_jn1.md)
  : CDC SARS-CoV-2 variant proportions: JN.1 emergence (US, 2023-2024)
- [`influenza_h3n2`](https://CuiweiG.github.io/lineagefreq/reference/influenza_h3n2.md)
  : Simulated influenza A/H3N2 clade frequency data
