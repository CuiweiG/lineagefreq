## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

## Test environments

* local Windows 10, R 4.5.3
* GitHub Actions: ubuntu-latest (R release, devel, oldrel-1)
* GitHub Actions: macos-latest (R release)
* GitHub Actions: windows-latest (R release)

## Package size

Installed size: ~1 MB.

## Reverse dependencies

New package, no reverse dependencies.

## Notes

* `Additional_repositories` field lists the Stan R-universe
  (`https://mc-stan.org/r-packages/`) because the optional
  Suggested dependency `cmdstanr` is not on CRAN. All
  CmdStan-dependent functionality is gated behind
  `lfq_stan_available()` and wrapped in `\donttest{}`
  or `skip_if_not()`.
