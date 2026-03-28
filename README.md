# lineagefreq

<!-- badges: start -->
<!-- badges: end -->

**Lineage Frequency Dynamics and Growth-Advantage Estimation from
Genomic Surveillance Counts**

lineagefreq models pathogen lineage frequency dynamics from genomic
surveillance count data. It provides a unified interface for
multiple modeling engines, rolling-origin backtesting, standardized
forecast scoring, emergence detection, and sequencing power analysis.

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("CuiweiG/lineagefreq")
```

## Quick example

```r
library(lineagefreq)

# Load built-in simulated SARS-CoV-2 data
data(sarscov2_us_2022)

# Prepare data
x <- lfq_data(sarscov2_us_2022,
              lineage = variant, date = date,
              count = count, total = total)

# Fit multinomial logistic regression
fit <- fit_model(x, engine = "mlr")

# Extract growth advantages
growth_advantage(fit, type = "relative_Rt", generation_time = 5)

# Forecast 4 weeks ahead
fc <- forecast(fit, horizon = 28)
autoplot(fc)

# Detect emerging lineages
summarize_emerging(x)
```

## Results on Real CDC Data

All figures below use **real surveillance data** from the U.S. CDC
([data.cdc.gov](https://data.cdc.gov/Laboratory-Surveillance/SARS-CoV-2-Variant-Proportions/jr58-6ysp),
public domain). Two independent epidemic waves are analyzed to
demonstrate robustness.

### Case Study 1 — JN.1 Emergence (2023–2024)

JN.1 rose from <1% to >80% in five months, displacing XBB-derived
lineages. MLR captures the complete sigmoid replacement curve.

![](man/figures/jn1_dynamics.png)

### Case Study 2 — BA.1 → BA.2 Replacement (2022)

The best-documented variant replacement event: four sequential
Omicron subvariant sweeps over six months.

![](man/figures/ba2_dynamics.png)

### Growth Advantage Estimation

Relative Rt estimates match published values: BA.2 = 1.34× vs BA.1
([Lyngse et al. 2022](https://doi.org/10.1038/s41467-022-33498-0),
published 1.3–1.5×); KP.3 = 1.25× vs JN.1.

![](man/figures/growth_advantage_comparison.png)

### Short-Term Frequency Forecast

Six-week projection with 95% prediction intervals from parametric
bootstrap (1,000 MVN draws from Fisher information matrix).

![](man/figures/forecast_plot.png)

### Rolling-Origin Forecast Accuracy

Out-of-sample evaluation: 3.8% MAE at 2-week horizon, 7.9% at
4-week horizon on the BA.2 period.

![](man/figures/backtest_plot.png)

## Features

- **Unified modeling**: `fit_model()` with engines `"mlr"`,
  `"hier_mlr"`, `"piantham"`.
- **Growth advantage**: four output types (growth rate, relative Rt,
  selection coefficient, doubling time).
- **Forecasting**: parametric simulation with prediction intervals.
- **Backtesting**: `backtest()` with rolling origins across multiple
  engines and horizons.
- **Scoring**: MAE, RMSE, coverage, weighted interval score via
  `score_forecasts()`.
- **Surveillance tools**: `summarize_emerging()`,
  `sequencing_power()`, `collapse_lineages()`.
- **broom compatible**: `tidy()`, `glance()`, `augment()`.
- **Publication-ready plots**: `autoplot()` for fits, forecasts,
  and backtest results.

## Supported pathogens

lineagefreq works with any pathogen that has variant/lineage-resolved
sequencing count data, including SARS-CoV-2, influenza, RSV, and
others.

## Citation

If you use lineagefreq in published work, please cite:

> [Package paper citation — to be updated upon publication]

## License

MIT
