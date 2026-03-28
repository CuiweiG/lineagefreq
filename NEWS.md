# lineagefreq 0.2.0 (development)

## New features

* **Pipe-friendly API**: `lfq_fit()`, `lfq_advantage()`, `lfq_forecast()`,
  `lfq_score()` enable tidyverse-style chaining with `|>`.
* **Engine registry**: `register_engine()` / `unregister_engine()` allow
  third-party packages to register custom modeling engines, similar to
  the parsnip engine system.
* **`lfq_summary()`**: One-row-per-lineage overview combining growth rates,
  confidence intervals, and relative Rt in a single tibble.
* **`as.data.frame.lfq_data()`**: Clean tibble export for interoperability.

## Improvements

* `fit_model()` now accepts both built-in and registered engine names.
* `lfq_engines()` lists all available engines including custom registrations.
* 251 unit tests (up from 243 in v0.1.0).

# lineagefreq 0.1.0

* Initial CRAN release.
