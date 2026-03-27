## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local Windows, R 4.5.3
* GitHub Actions: ubuntu-latest (R release)

## Package size

Installed size: ~1 MB (well under the 5 MB limit).

## Reverse dependencies

New package, no reverse dependencies.

## Notes

* `Additional_repositories` lists the Stan R-universe for the
  optional `cmdstanr` Suggested dependency (not on CRAN).
  Bayesian engine functionality is gated behind
  `lfq_stan_available()` and all related examples/tests are
  wrapped in `\donttest{}` / `skip_if_not()`.
