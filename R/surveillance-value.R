#' Expected Value of Information for genomic surveillance
#'
#' Quantifies the marginal decision value of sequencing additional
#' samples from each region or stratum, based on the current
#' posterior uncertainty of variant frequency estimates. The EVOI
#' captures how much the expected estimation error decreases per
#' additional sequence, enabling cost-effective resource allocation.
#'
#' @param fit An \code{lfq_fit} object from \code{\link{fit_model}}.
#' @param n_current Integer vector of current sample sizes per
#'   stratum. If a single integer, assumed equal across all lineages.
#' @param n_additional Integer vector of candidate additional sample
#'   sizes to evaluate. Default \code{seq(10, 200, by = 10)}.
#' @param target_lineage Optional character; compute EVOI
#'   specifically for this lineage. Default \code{NULL} evaluates
#'   across all non-pivot lineages.
#'
#' @return An \code{evoi} S3 class with components:
#'   \describe{
#'     \item{values}{Tibble with \code{n_additional}, \code{evoi},
#'       \code{marginal_evoi} columns.}
#'     \item{current_uncertainty}{Numeric; current estimation
#'       variance (sum of growth rate SEs squared).}
#'     \item{target_lineage}{Character or NULL.}
#'   }
#'
#' @details
#' The EVOI is computed as the expected reduction in mean squared
#' estimation error for lineage frequency when \eqn{n} additional
#' sequences are observed. Under the multinomial likelihood with
#' Gaussian approximation to the posterior, the variance of the
#' frequency estimate scales as \eqn{p(1-p) / n}, and additional
#' samples reduce this by the factor \eqn{n_0 / (n_0 + n_{\text{add}})}.
#'
#' The marginal EVOI (the value of one additional sequence) is
#' the derivative of the EVOI curve. It decreases monotonically
#' with sample size, exhibiting the diminishing returns
#' characteristic of information-theoretic quantities.
#'
#' @seealso \code{\link{adaptive_design}} for dynamic allocation,
#'   \code{\link{sequencing_power}} for static sample size planning.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.3, "B" = 0.9),
#'   n_timepoints = 15, seed = 1)
#' fit <- fit_model(sim, engine = "mlr")
#' ev <- surveillance_value(fit, n_current = 500)
#' ev
#' }
#'
#' @export
surveillance_value <- function(fit,
                               n_current    = 100L,
                               n_additional = seq(10L, 200L, by = 10L),
                               target_lineage = NULL) {

  if (!inherits(fit, "lfq_fit")) {
    cli::cli_abort("{.arg fit} must be an {.cls lfq_fit} object.")
  }

  lineages  <- fit$lineages
  pivot     <- fit$pivot
  non_pivot <- setdiff(lineages, pivot)
  vcov_mat  <- fit$vcov_matrix

  if (!is.null(target_lineage)) {
    if (!target_lineage %in% non_pivot) {
      cli::cli_abort(
        "{.arg target_lineage} must be a non-pivot lineage."
      )
    }
    non_pivot <- target_lineage
  }

  # Current uncertainty: sum of growth rate variances
  delta_names <- paste0("delta_", non_pivot)
  current_var <- sum(vapply(delta_names, function(nm) {
    if (nm %in% rownames(vcov_mat)) vcov_mat[nm, nm] else 0
  }, numeric(1L)))

  n0 <- as.numeric(n_current[1])

  # EVOI: expected variance reduction
  evoi_vals <- vapply(n_additional, function(n_add) {
    # Variance scales as 1/n; adding n_add reduces by factor n0/(n0+n_add)
    reduction_factor <- n0 / (n0 + n_add)
    current_var * (1 - reduction_factor)
  }, numeric(1L))

  # Marginal EVOI: discrete derivative
  marginal <- c(evoi_vals[1],
                diff(evoi_vals) / diff(as.numeric(n_additional)))

  values <- tibble::tibble(
    n_additional  = as.integer(n_additional),
    evoi          = evoi_vals,
    marginal_evoi = marginal
  )

  structure(
    list(
      values              = values,
      current_uncertainty = current_var,
      target_lineage      = target_lineage
    ),
    class = "evoi"
  )
}


#' @export
print.evoi <- function(x, ...) {
  cli::cli_h3("Expected Value of Information")
  cli::cli_text(
    "Current estimation variance: {round(x$current_uncertainty, 6)}"
  )
  if (!is.null(x$target_lineage)) {
    cli::cli_text("Target lineage: {x$target_lineage}")
  }
  n50 <- x$values$n_additional[which.min(abs(x$values$evoi -
    0.5 * max(x$values$evoi)))]
  cli::cli_text(
    "50% variance reduction at ~{n50} additional sequences"
  )
  invisible(x)
}


#' Plot EVOI curve
#' @param x An \code{evoi} object.
#' @param ... Unused.
#' @return A ggplot object.
#' @export
plot.evoi <- function(x, ...) {
  ggplot2::ggplot(x$values,
    ggplot2::aes(x = .data$n_additional, y = .data$evoi)) +
    ggplot2::geom_line(colour = "#2166ac", linewidth = 0.9) +
    ggplot2::geom_point(colour = "#2166ac", size = 1.5) +
    ggplot2::labs(
      x = "Additional sequences",
      y = "Expected variance reduction",
      title = "Diminishing returns of additional sequencing"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}
