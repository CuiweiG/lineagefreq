# Changelog

## lineagefreq 0.2.0

CRAN release: 2026-04-03

### New features

- **Pipe-friendly API**:
  [`lfq_fit()`](https://cuiweig.github.io/lineagefreq/reference/lfq_fit.md),
  [`lfq_advantage()`](https://cuiweig.github.io/lineagefreq/reference/lfq_advantage.md),
  [`lfq_forecast()`](https://cuiweig.github.io/lineagefreq/reference/lfq_forecast.md),
  [`lfq_score()`](https://cuiweig.github.io/lineagefreq/reference/lfq_score.md)
  enable tidyverse-style chaining with `|>`.
- **Engine registry**:
  [`register_engine()`](https://cuiweig.github.io/lineagefreq/reference/register_engine.md)
  /
  [`unregister_engine()`](https://cuiweig.github.io/lineagefreq/reference/unregister_engine.md)
  allow third-party packages to register custom modeling engines,
  similar to the parsnip engine system.
- **[`lfq_summary()`](https://cuiweig.github.io/lineagefreq/reference/lfq_summary.md)**:
  One-row-per-lineage overview combining growth rates, confidence
  intervals, and relative Rt in a single tibble.
- **[`as.data.frame.lfq_data()`](https://cuiweig.github.io/lineagefreq/reference/as.data.frame.lfq_data.md)**:
  Clean tibble export for interoperability.

### Improvements

- [`fit_model()`](https://cuiweig.github.io/lineagefreq/reference/fit_model.md)
  now accepts both built-in and registered engine names.
- [`lfq_engines()`](https://cuiweig.github.io/lineagefreq/reference/lfq_engines.md)
  lists all available engines including custom registrations.
- 251 unit tests (up from 243 in v0.1.0).

## lineagefreq 0.1.0

- Initial CRAN release.
