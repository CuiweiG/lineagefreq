skip_if_no_stan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not(
    lineagefreq::lfq_stan_available(),
    message = "CmdStan not available"
  )
}
