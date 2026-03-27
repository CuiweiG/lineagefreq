## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

## Test environments

* local Windows 10, R 4.5.3
* GitHub Actions: ubuntu-latest (R release, R-devel, R oldrel-1)
* GitHub Actions: macos-latest (R release)
* GitHub Actions: windows-latest (R release)

## Package size

Installed size: approximately 1 MB (tarball: 103 KB).

## Reverse dependencies

This is a new submission. No reverse dependencies.

## Additional notes

* The 'Additional_repositories' field lists the Stan R-universe
  (https://mc-stan.org/r-packages/) for the optional Suggested
  dependency 'cmdstanr', which is not on CRAN. All functionality
  requiring 'CmdStan' is gated behind lfq_stan_available() checks
  and corresponding examples/tests use \donttest{} and
  skip_if_not() respectively.
