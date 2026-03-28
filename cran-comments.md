## Resubmission

This is a resubmission. Changes from previous submission:

* Removed non-standard top-level file (submit.R).
* Words flagged as possibly misspelled in DESCRIPTION are correct:
  Abousamra, Figgins, Bedford (author surnames); Piantham (method
  name from published literature); backtesting (standard term in
  forecast evaluation).
* Shortened Title to 59 characters.

## Package motivation

No CRAN package currently provides a unified interface for modeling
pathogen variant/lineage frequency dynamics from genomic surveillance
count data. Existing tools (e.g., Nextstrain's evofr) are not on
CRAN and lack built-in forecast evaluation. lineagefreq fills this
gap by providing:

* Multinomial logistic regression with multiple engine backends
* Growth advantage estimation validated against published values
  (Lyngse et al. 2022, Abousamra et al. 2024)
* Probabilistic forecasting with prediction intervals
* Rolling-origin backtesting framework for honest accuracy assessment

The package targets public health genomic surveillance teams working
with SARS-CoV-2, influenza, RSV, and other variant-resolved pathogens.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

## Test environments

* local Windows 10, R 4.5.3
* GitHub Actions: ubuntu-latest (R release, R oldrel-1)
* GitHub Actions: macOS-latest (R release)
* GitHub Actions: windows-latest (R release)

## Reverse dependencies

New package, no reverse dependencies.
