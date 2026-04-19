# lineagefreq: Lineage Frequency Dynamics from Genomic Surveillance Counts

Models pathogen lineage frequency dynamics from genomic surveillance
count data. Provides a unified
[`fit_model()`](https://cuiweig.github.io/lineagefreq/reference/fit_model.md)
interface with multiple engines (MLR, hierarchical MLR, Piantham),
rolling-origin backtesting via
[`backtest()`](https://cuiweig.github.io/lineagefreq/reference/backtest.md),
standardized scoring via
[`score_forecasts()`](https://cuiweig.github.io/lineagefreq/reference/score_forecasts.md),
emergence detection, and sequencing power analysis.

## Quick start

    x <- lfq_data(my_counts, lineage = lineage, date = date, count = n)
    fit <- fit_model(x, engine = "mlr")
    growth_advantage(fit, generation_time = 5)

## See also

Useful links:

- <https://github.com/CuiweiG/lineagefreq>

- <https://CuiweiG.github.io/lineagefreq>

- Report bugs at <https://github.com/CuiweiG/lineagefreq/issues>

## Author

**Maintainer**: Cuiwei Gao <48gaocuiwei@gmail.com> \[copyright holder\]
