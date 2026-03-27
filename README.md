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
