#' Detection horizon for an emerging variant
#'
#' Given current sequencing capacity and a hypothetical variant at
#' initial prevalence \eqn{p_0} growing at rate \eqn{\delta},
#' estimates the number of surveillance periods (weeks) until the
#' probability of detection exceeds a specified threshold.
#'
#' @param initial_prev Initial prevalence of the emerging variant.
#'   Default 0.001 (0.1 percent).
#' @param growth_rate Per-week multiplicative growth rate on the
#'   frequency scale. Default 1.3 (30 percent per week).
#' @param n_per_period Sequences collected per surveillance period.
#'   Default 500.
#' @param n_periods Maximum number of periods to evaluate. Default
#'   26 (6 months of weekly data).
#' @param detection_threshold Minimum number of sequences of the
#'   target lineage required to declare detection. Default 1.
#' @param confidence Detection probability threshold. Default 0.95.
#'
#' @return A tibble with columns \code{period}, \code{prevalence},
#'   \code{detection_prob} (cumulative detection probability), and
#'   \code{detected} (logical). An attribute \code{weeks_to_detection}
#'   contains the first period where detection probability exceeds
#'   the confidence threshold, or \code{NA} if not reached.
#'
#' @details
#' At each period \eqn{t}, the variant prevalence is modelled as
#' logistic growth: \eqn{p(t) = p_0 \cdot r^t / (1 - p_0 +
#' p_0 \cdot r^t)} where \eqn{r} is the per-period growth rate.
#' The probability of detecting at least \eqn{k} sequences in
#' \eqn{n} draws at prevalence \eqn{p} is
#' \eqn{1 - F_{\text{Binom}}(k-1; n, p)}.
#' Cumulative detection is \eqn{1 - \prod_{\tau=1}^{t}(1 - P_\tau)}.
#'
#' @seealso \code{\link{sequencing_power}} for static sample size
#'   calculations, \code{\link{alert_threshold}} for sequential
#'   detection.
#'
#' @examples
#' dh <- detection_horizon(initial_prev = 0.005, growth_rate = 1.2,
#'   n_per_period = 300)
#' dh
#'
#' @export
detection_horizon <- function(initial_prev        = 0.001,
                              growth_rate          = 1.3,
                              n_per_period         = 500L,
                              n_periods            = 26L,
                              detection_threshold  = 1L,
                              confidence           = 0.95) {

  if (initial_prev <= 0 || initial_prev >= 1) {
    cli::cli_abort("{.arg initial_prev} must be in (0, 1).")
  }
  if (growth_rate <= 0) {
    cli::cli_abort("{.arg growth_rate} must be positive.")
  }

  # Logistic growth: p(t) = p0*r^t / (1 - p0 + p0*r^t)
  periods <- seq_len(n_periods)
  prevalence <- vapply(periods, function(t) {
    numerator <- initial_prev * growth_rate^t
    numerator / (1 - initial_prev + numerator)
  }, numeric(1L))

  # Per-period detection probability
  per_period_detect <- stats::pbinom(
    detection_threshold - 1L, size = n_per_period, prob = prevalence,
    lower.tail = FALSE
  )

  # Cumulative detection: 1 - prod(1 - P_t)
  cum_miss <- cumprod(1 - per_period_detect)
  cum_detect <- 1 - cum_miss

  detected <- cum_detect >= confidence
  weeks_to_detect <- if (any(detected)) min(which(detected)) else NA_integer_

  result <- tibble::tibble(
    period         = periods,
    prevalence     = prevalence,
    detection_prob = cum_detect,
    detected       = detected
  )

  attr(result, "weeks_to_detection") <- weeks_to_detect
  result
}
