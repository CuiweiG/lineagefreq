#' Population-level selective pressure from variant dynamics
#'
#' Computes a population-level selective pressure metric from
#' genomic surveillance data alone, without requiring case counts
#' or epidemiological data. The metric quantifies how rapidly the
#' variant landscape is shifting and serves as an early warning
#' signal for epidemic growth that is robust to case underreporting.
#'
#' The approach follows the framework of Figgins and Bedford (2025),
#' where selective pressure is defined as the variance-weighted mean
#' growth rate across circulating lineages:
#' \deqn{S(t) = \sum_v p_v(t) \cdot \delta_v}
#' This represents the expected rate at which the population-average
#' fitness is increasing, measured entirely from sequence data.
#'
#' @param fit An \code{lfq_fit} object from \code{\link{fit_model}}.
#' @param method Aggregation method: \code{"mean"} (default) for
#'   the frequency-weighted mean growth rate, or \code{"variance"}
#'   for the frequency-weighted variance of growth rates (higher
#'   variance indicates stronger selection).
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{date}{Date.}
#'     \item{pressure}{Selective pressure value.}
#'     \item{dominant_lineage}{Lineage with highest frequency.}
#'     \item{dominant_freq}{Frequency of dominant lineage.}
#'   }
#'
#' @details
#' When \code{method = "mean"}, the metric is positive when fitter-
#' than-average lineages are increasing in frequency, indicating
#' population-level adaptation. A sustained positive value precedes
#' epidemic growth because it means the effective reproduction number
#' of the average circulating virus is increasing.
#'
#' When \code{method = "variance"}, the metric captures the
#' heterogeneity of fitness across co-circulating lineages. High
#' variance indicates strong directional selection; low variance
#' indicates near-neutral drift.
#'
#' This metric requires only genomic surveillance data. It does not
#' require case counts, hospitalisations, or wastewater data, making
#' it applicable in settings where epidemiological reporting is
#' incomplete or delayed.
#'
#' @references
#' Figgins MD, Bedford T (2025). Jointly modeling variant frequencies
#' and case counts to estimate relative variant severity.
#' \emph{medRxiv}. \doi{10.1101/2024.12.02.24318334}
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 4,
#'   advantages = c("A" = 1.4, "B" = 1.1, "C" = 0.8),
#'   n_timepoints = 15, seed = 1)
#' fit <- fit_model(sim, engine = "mlr")
#' sp <- selective_pressure(fit)
#' sp
#' }
#'
#' @export
selective_pressure <- function(fit,
                               method = c("mean", "variance")) {

  if (!inherits(fit, "lfq_fit")) {
    cli::cli_abort("{.arg fit} must be an {.cls lfq_fit} object.")
  }

  method    <- match.arg(method)
  lineages  <- fit$lineages
  deltas    <- fit$growth_rates
  fitted_df <- fit$fitted_values
  dates     <- sort(unique(fitted_df$.date))

  rows <- list()

  for (d in dates) {
    d <- as.Date(d, origin = "1970-01-01")
    sub <- fitted_df[fitted_df$.date == d, ]

    freqs <- stats::setNames(sub$.fitted_freq, sub$.lineage)
    freqs <- freqs[lineages]

    if (method == "mean") {
      # Frequency-weighted mean growth rate
      pressure <- sum(freqs * deltas[lineages], na.rm = TRUE)
    } else {
      # Frequency-weighted variance
      mean_delta <- sum(freqs * deltas[lineages], na.rm = TRUE)
      pressure <- sum(freqs * (deltas[lineages] - mean_delta)^2,
                      na.rm = TRUE)
    }

    dominant <- names(which.max(freqs))
    dominant_freq <- unname(freqs[dominant])

    rows <- c(rows, list(tibble::tibble(
      date           = d,
      pressure       = pressure,
      dominant_lineage = dominant,
      dominant_freq  = dominant_freq
    )))
  }

  dplyr::bind_rows(rows)
}
