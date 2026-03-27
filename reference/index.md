# Package index

## Data preparation

- [`lfq_data()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
  : Create a lineage frequency data object
- [`is_lfq_data()`](https://CuiweiG.github.io/lineagefreq/reference/is_lfq_data.md)
  : Test if an object is an lfq_data object
- [`collapse_lineages()`](https://CuiweiG.github.io/lineagefreq/reference/collapse_lineages.md)
  : Collapse rare lineages into an aggregate group
- [`filter_sparse()`](https://CuiweiG.github.io/lineagefreq/reference/filter_sparse.md)
  : Filter sparse time points and lineages
- [`simulate_dynamics()`](https://CuiweiG.github.io/lineagefreq/reference/simulate_dynamics.md)
  : Simulate lineage frequency dynamics

## Model fitting

- [`fit_model()`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md)
  : Fit a lineage frequency model

## Inference and forecasting

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

## Backtesting and evaluation

- [`backtest()`](https://CuiweiG.github.io/lineagefreq/reference/backtest.md)
  : Rolling-origin backtesting of lineage frequency models
- [`score_forecasts()`](https://CuiweiG.github.io/lineagefreq/reference/score_forecasts.md)
  : Score backtest forecast accuracy
- [`compare_models()`](https://CuiweiG.github.io/lineagefreq/reference/compare_models.md)
  : Compare model engines from backtest scores
- [`plot_backtest()`](https://CuiweiG.github.io/lineagefreq/reference/plot_backtest.md)
  : Plot backtest scores

## S3 methods

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
- [`autoplot(`*`<lfq_fit>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/autoplot.lfq_fit.md)
  : Plot lineage frequency model results
- [`autoplot(`*`<lfq_forecast>`*`)`](https://CuiweiG.github.io/lineagefreq/reference/autoplot.lfq_forecast.md)
  : Plot a lineage frequency forecast

## Datasets

- [`sarscov2_us_2022`](https://CuiweiG.github.io/lineagefreq/reference/sarscov2_us_2022.md)
  : Simulated SARS-CoV-2 variant frequency data (US, 2022)
- [`influenza_h3n2`](https://CuiweiG.github.io/lineagefreq/reference/influenza_h3n2.md)
  : Simulated influenza A/H3N2 clade frequency data

## Utilities

- [`lfq_stan_available()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_stan_available.md)
  : Check if CmdStan backend is available
