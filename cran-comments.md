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

Installed size: approximately 1 MB (tarball: ~103 KB).

## Reverse dependencies

This is a new submission. No reverse dependencies.

## Additional notes

* Five modeling engines are provided: 'mlr', 'hier_mlr', 'piantham'
  (frequentist, no external dependencies), and 'fga'/'garw'
  (Bayesian, require 'CmdStan').

* The 'Additional_repositories' field lists the Stan R-universe
  (https://mc-stan.org/r-packages/) for the optional Suggested

  dependency 'cmdstanr', which is not on CRAN. All functionality
  requiring 'CmdStan' is gated behind `lfq_stan_available()` and
  corresponding examples use `\donttest{}`, tests use `skip_if_not()`.

* 197 unit tests pass, 6 are skipped (require 'CmdStan').
  Test coverage: 91%.
