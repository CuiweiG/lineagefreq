#' Sequential detection of emerging variants
#'
#' Applies a sequential probability ratio test (SPRT) or CUSUM
#' procedure to lineage frequency data, determining when accumulated
#' evidence is sufficient to declare a variant "emerging" rather
#' than sampling noise. Controls the false alarm rate while
#' minimising detection delay.
#'
#' @param data An \code{lfq_data} object.
#' @param method Detection method: \code{"sprt"} (default) for
#'   the sequential probability ratio test, or \code{"cusum"} for
#'   the cumulative sum control chart.
#' @param alpha False alarm probability. Default 0.05.
#' @param beta Missed detection probability. Default 0.10.
#' @param delta_0 Null hypothesis growth rate (no emergence).
#'   Default 0 (frequency is stable).
#' @param delta_1 Alternative hypothesis growth rate (emergence).
#'   Default 0.03 (3 percent per-week increase on logit scale).
#' @param threshold CUSUM decision threshold. Default 5.0.
#'   Only used when \code{method = "cusum"}.
#'
#' @return A tibble with columns \code{lineage}, \code{date},
#'   \code{statistic} (log-likelihood ratio or CUSUM value),
#'   \code{alert} (logical), \code{direction} (emerging/declining/
#'   stable), and \code{confidence} (1 - alpha).
#'
#' @details
#' \strong{SPRT} (Wald, 1945) computes the log-likelihood ratio
#' between the alternative (lineage is growing at rate
#' \eqn{\delta_1}) and the null (frequency is stable). The test
#' stops when the cumulative log-ratio crosses the upper boundary
#' \eqn{B = \log((1-\beta)/\alpha)} (declare emerging) or the
#' lower boundary \eqn{A = \log(\beta/(1-\alpha))} (declare stable).
#'
#' \strong{CUSUM} accumulates deviations from expected frequency
#' under the null: \eqn{S_t = \max(0, S_{t-1} + (x_t - k))}
#' where \eqn{x_t} is the observed frequency change and \eqn{k}
#' is the allowance (half the shift to detect). An alert is raised
#' when \eqn{S_t > h}.
#'
#' @references
#' Wald A (1945). Sequential tests of statistical hypotheses.
#' \emph{Annals of Mathematical Statistics}, 16(2), 117--186.
#'
#' @seealso \code{\link{summarize_emerging}} for non-sequential
#'   trend tests, \code{\link{detection_horizon}} for prospective
#'   power analysis.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.5, "B" = 0.8),
#'   n_timepoints = 15, seed = 1)
#' alerts <- alert_threshold(sim)
#' alerts
#' }
#'
#' @export
alert_threshold <- function(data,
                            method    = c("sprt", "cusum"),
                            alpha     = 0.05,
                            beta      = 0.10,
                            delta_0   = 0,
                            delta_1   = 0.03,
                            threshold = 5.0) {

  if (!is_lfq_data(data)) {
    cli::cli_abort("{.arg data} must be an {.cls lfq_data} object.")
  }

  method   <- match.arg(method)
  lineages <- attr(data, "lineages")
  dates    <- sort(unique(data$.date))

  results <- list()

  for (lin in lineages) {
    lin_data <- data[data$.lineage == lin, ]
    lin_data <- lin_data[order(lin_data$.date), ]

    if (nrow(lin_data) < 3L) {
      results <- c(results, list(tibble::tibble(
        lineage    = lin,
        date       = max(lin_data$.date),
        statistic  = 0,
        alert      = FALSE,
        direction  = "stable",
        confidence = 1 - alpha
      )))
      next
    }

    freqs <- lin_data$.freq

    if (method == "sprt") {
      # SPRT boundaries
      upper_b <- log((1 - beta) / alpha)
      lower_b <- log(beta / (1 - alpha))

      # Compute incremental log-likelihood ratios
      cum_llr <- 0
      alert_raised <- FALSE
      alert_date   <- max(lin_data$.date)

      for (i in 2:length(freqs)) {
        # Change in logit-frequency
        logit_prev <- stats::qlogis(pmin(pmax(freqs[i - 1], 0.001), 0.999))
        logit_curr <- stats::qlogis(pmin(pmax(freqs[i], 0.001), 0.999))
        x <- logit_curr - logit_prev

        # Log-likelihood ratio: N(delta_1, 1) vs N(delta_0, 1)
        llr <- delta_1 * x - 0.5 * delta_1^2 -
               (delta_0 * x - 0.5 * delta_0^2)
        cum_llr <- cum_llr + llr

        if (cum_llr >= upper_b) {
          alert_raised <- TRUE
          alert_date   <- lin_data$.date[i]
          break
        }
        if (cum_llr <= lower_b) {
          break
        }
      }

      direction <- if (alert_raised) {
        "emerging"
      } else if (cum_llr <= lower_b) {
        "stable"
      } else {
        "inconclusive"
      }

      results <- c(results, list(tibble::tibble(
        lineage    = lin,
        date       = alert_date,
        statistic  = cum_llr,
        alert      = alert_raised,
        direction  = direction,
        confidence = 1 - alpha
      )))

    } else {
      # CUSUM
      k <- delta_1 / 2  # allowance
      cusum <- 0
      alert_raised <- FALSE
      alert_date   <- max(lin_data$.date)

      for (i in 2:length(freqs)) {
        logit_prev <- stats::qlogis(pmin(pmax(freqs[i - 1], 0.001), 0.999))
        logit_curr <- stats::qlogis(pmin(pmax(freqs[i], 0.001), 0.999))
        x <- logit_curr - logit_prev

        cusum <- max(0, cusum + x - k)

        if (cusum >= threshold) {
          alert_raised <- TRUE
          alert_date   <- lin_data$.date[i]
          break
        }
      }

      results <- c(results, list(tibble::tibble(
        lineage    = lin,
        date       = alert_date,
        statistic  = cusum,
        alert      = alert_raised,
        direction  = if (alert_raised) "emerging" else "stable",
        confidence = 1 - alpha
      )))
    }
  }

  dplyr::bind_rows(results)
}
