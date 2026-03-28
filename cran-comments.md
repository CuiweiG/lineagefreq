## Resubmission

This is a resubmission. Changes from the previous submission:

* Removed `submit.R` from top-level directory.
* The words flagged as possibly misspelled in DESCRIPTION are all
  correct: Abousamra, Figgins, and Bedford are author surnames from
  the cited reference; Piantham is a method name from published
  literature; backtesting is a standard term in forecast evaluation.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

  Suggests or Enhances not in mainstream repositories:
    cmdstanr
  Availability using Additional_repositories specification:
    cmdstanr   yes   https://mc-stan.org/r-packages/

## Test environments

* local Windows 10, R 4.5.3
* GitHub Actions: ubuntu-latest (R release, R oldrel-1)
* GitHub Actions: macOS-latest (R release)
* GitHub Actions: windows-latest (R release)

## Package size

Installed size: approximately 1 MB (tarball: ~1.4 MB including
vignettes and PNG figures).

## Reverse dependencies

This is a new submission. No reverse dependencies.

## Additional notes

* `cmdstanr` is listed in Suggests with `Additional_repositories`
  pointing to `https://mc-stan.org/r-packages/`. It is only checked
  at runtime via `lfq_stan_available()`. No tests, examples, or
  vignettes depend on it.
