## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

## Test environments

* local Windows 10, R 4.5.3
* GitHub Actions: ubuntu-latest (R release, devel, oldrel-1)
* GitHub Actions: macos-latest (R release)
* GitHub Actions: windows-latest (R release)

## Package size

Installed size: ~1 MB (tarball: 103 KB).

## Reverse dependencies

New package, no reverse dependencies.

## Notes

* `cmdstanr` is listed in Suggests with `Additional_repositories`
  pointing to `https://mc-stan.org/r-packages/`. It is only used
  in `lfq_stan_available()` to check backend availability for
  future Bayesian engines (v0.2). No tests or examples depend on it.
