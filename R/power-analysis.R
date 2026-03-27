#' Sequencing power analysis
#'
#' Estimates the minimum number of sequences needed to detect a
#' lineage at a given frequency with specified precision.
#'
#' @param target_precision Desired half-width of the frequency
#'   confidence interval. Default 0.05 (plus/minus 5 percentage points).
#' @param current_freq True or assumed frequency of the target
#'   lineage. Can be a vector for multiple scenarios.
#'   Default 0.02 (2%).
#' @param ci_level Confidence level. Default 0.95.
#'
#' @return A tibble with columns: `current_freq`,
#'   `target_precision`, `required_n`, `ci_level`.
#'
#' @details
#' Uses the normal approximation to the binomial:
#' \deqn{n = z^2 \cdot p(1-p) / E^2}
#' where z is the critical value, p is frequency, E is precision.
#'
#' @examples
#' # How many sequences to estimate a 2% lineage within +/-5%?
#' sequencing_power()
#'
#' # Multiple scenarios
#' sequencing_power(current_freq = c(0.01, 0.02, 0.05, 0.10))
#'
#' @export
sequencing_power <- function(target_precision = 0.05,
                             current_freq     = 0.02,
                             ci_level         = 0.95) {

  if (any(target_precision <= 0)) {
    cli::cli_abort("{.arg target_precision} must be positive.")
  }
  if (any(current_freq <= 0 | current_freq >= 1)) {
    cli::cli_abort("{.arg current_freq} must be in (0, 1).")
  }
  assert_prob(ci_level, "ci_level")

  z <- stats::qnorm((1 + ci_level) / 2)

  tibble::tibble(
    current_freq     = current_freq,
    target_precision = target_precision,
    required_n       = ceiling(z^2 * current_freq * (1 - current_freq) /
                                 target_precision^2),
    ci_level         = ci_level
  )
}
