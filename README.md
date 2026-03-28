# lineagefreq

*Lineage Frequency Dynamics and Growth-Advantage Estimation from Genomic Surveillance Counts*

<!-- badges: start -->
[![R-CMD-check](https://github.com/CuiweiG/lineagefreq/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CuiweiG/lineagefreq/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://img.shields.io/badge/CRAN-submitted-orange.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R ≥ 4.1.0](https://img.shields.io/badge/R-%E2%89%A5%204.1.0-brightgreen.svg)](https://cran.r-project.org/)
<!-- badges: end -->

An R package for modeling pathogen lineage frequencies, estimating
growth advantages, and forecasting variant replacement dynamics from
genomic surveillance counts.

## Installation

```r
# install.packages("pak")
pak::pak("CuiweiG/lineagefreq")

# Or with devtools:
# devtools::install_github("CuiweiG/lineagefreq")
```

## Quick example

```r
library(lineagefreq)
library(ggplot2)

data(cdc_sarscov2_jn1)
x <- lfq_data(cdc_sarscov2_jn1,
              lineage = lineage, date = date, count = count)

fit <- fit_model(x, engine = "mlr")
growth_advantage(fit, type = "relative_Rt", generation_time = 5)

fc <- forecast(fit, horizon = 28)
autoplot(fc)
```

## Real-Data Case Studies

Figures below use **real U.S. CDC surveillance data**
([data.cdc.gov/jr58-6ysp](https://data.cdc.gov/Laboratory-Surveillance/SARS-CoV-2-Variant-Proportions/jr58-6ysp),
public domain). Two independent epidemic waves illustrate model
behavior across distinct replacement settings.

Data accessed 2026-03-28. Lineages below 5% peak frequency collapsed
to "Other." Reproducible scripts: `data-raw/prepare_cdc_data.R` and
`data-raw/prepare_ba2_data.R`.

### Variant Replacement Dynamics

**JN.1 emergence (Oct 2023 – Mar 2024):** MLR recovers the observed
replacement trajectory from <1% to >80%.

![](man/figures/jn1_dynamics.png)

**BA.1 → BA.2 period (Dec 2021 – Jun 2022):** A well-characterized
Omicron replacement wave with four sequential subvariant sweeps.

![](man/figures/ba2_dynamics.png)

### Growth Advantage Estimation

Relative Rt estimates are consistent with published values:
BA.2 = 1.34× vs BA.1
([Lyngse et al. 2022](https://doi.org/10.1038/s41467-022-33498-0),
published 1.3–1.5×); KP.3 = 1.36× vs JN.1. Generation times:
3.2 days for Omicron BA.* subvariants
([Du et al. 2022](https://doi.org/10.3201/eid2806.220158));
5.0 days for JN/KP lineages.

![](man/figures/growth_advantage_comparison.png)

### Frequency Forecast

Six-week projection with 95% marginal prediction intervals
(pointwise, not simultaneous). Uncertainty reflects parameter
estimation error (MVN from Fisher information) and multinomial
sampling noise (n_eff = 100 sequences/period). See figure caption
for full methodological notes.

![](man/figures/forecast_plot.png)

### Forecast Accuracy

Rolling-origin out-of-sample evaluation on the BA.2 period:
3.9% MAE at 2-week horizon, 7.7% at 4-week horizon.

![](man/figures/backtest_plot.png)

## Features

**Model fitting**
- `fit_model()` with engines `"mlr"`, `"hier_mlr"`, `"piantham"`,
  `"fga"`, `"garw"` (Bayesian engines require
  ['CmdStan'](https://mc-stan.org/cmdstanr/))

**Inference**
- Growth advantage in four scales: growth rate, relative Rt,
  selection coefficient, doubling time

**Forecasting**
- Probabilistic frequency forecasts with parametric simulation
  and configurable sampling noise

**Evaluation**
- Rolling-origin backtesting via `backtest()` with standardized
  scoring (MAE, RMSE, coverage, WIS) via `score_forecasts()`

**Surveillance utilities**
- `summarize_emerging()`: binomial GLM trend tests per lineage
- `sequencing_power()`: minimum sample size for detection
- `collapse_lineages()`, `filter_sparse()`: preprocessing

**Visualization**
- `autoplot()` methods for fits, forecasts, and backtest summaries
- Publication-quality output with colorblind-safe palettes

**Interoperability**
- broom-compatible: `tidy()`, `glance()`, `augment()`
- `as_lfq_data()` generic for extensible data import
- `read_lineage_counts()` for CSV input

## Supported pathogens

Any pathogen with variant/lineage-resolved sequencing count data:
SARS-CoV-2, influenza, RSV, mpox, and others.

## Citation

```r
citation("lineagefreq")
```

A software paper and Zenodo DOI will be added upon publication.

## License

MIT
