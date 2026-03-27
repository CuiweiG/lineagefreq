#' @noRd
assert_col <- function(data, col, arg = deparse(substitute(data))) {
  if (!col %in% names(data)) {
    cli::cli_abort("Column {.field {col}} not found in {.arg {arg}}.",
                   call = rlang::caller_env())
  }
}

#' @noRd
assert_pos_int <- function(x, nm = deparse(substitute(x))) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      x < 1 || x != round(x)) {
    cli::cli_abort("{.arg {nm}} must be a positive integer.",
                   call = rlang::caller_env())
  }
}

#' @noRd
assert_prob <- function(x, nm = deparse(substitute(x))) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      x <= 0 || x >= 1) {
    cli::cli_abort("{.arg {nm}} must be in (0, 1).",
                   call = rlang::caller_env())
  }
}

#' @noRd
assert_date <- function(x, nm = deparse(substitute(x))) {
  if (!inherits(x, "Date")) {
    cli::cli_abort("{.arg {nm}} must be Date. Got {.cls {class(x)}}.",
                   call = rlang::caller_env())
  }
}

#' Check if CmdStan backend is available
#'
#' Returns `TRUE` if CmdStanR and CmdStan are installed.
#' Bayesian engines require this.
#'
#' @return Logical scalar.
#'
#' @examples
#' lfq_stan_available()
#'
#' @export
lfq_stan_available <- function() {
  requireNamespace("cmdstanr", quietly = TRUE) &&
    tryCatch(nchar(cmdstanr::cmdstan_path()) > 0,
             error = function(e) FALSE)
}

#' Log-sum-exp (numerically stable)
#' @noRd
log_sum_exp <- function(x) {
  mx <- max(x)
  mx + log(sum(exp(x - mx)))
}
