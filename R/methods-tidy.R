#' Tidy an lfq_fit object
#'
#' Converts model results to a tidy tibble, compatible with the
#' broom package ecosystem.
#'
#' @param x An [lfq_fit] object.
#' @param conf.int Include confidence intervals? Default `TRUE`.
#' @param conf.level Confidence level. Default 0.95.
#' @param ... Ignored.
#'
#' @return A tibble with columns: `lineage`, `term`, `estimate`,
#'   `std.error`, `conf.low`, `conf.high`.
#'
#' @examples
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
#' fit <- fit_model(sim)
#' tidy.lfq_fit(fit)
#'
#' @export
tidy.lfq_fit <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  non_pivot <- setdiff(x$lineages, x$pivot)
  vcov_mat  <- x$vcov_matrix
  z         <- stats::qnorm((1 + conf.level) / 2)

  rows <- list()
  for (v in non_pivot) {
    a_name <- paste0("alpha_", v)
    a_idx  <- match(a_name, colnames(vcov_mat))
    a_est  <- x$intercepts[v]
    a_se   <- if (!is.na(a_idx)) sqrt(max(vcov_mat[a_idx, a_idx], 0)) else NA_real_

    d_name <- paste0("delta_", v)
    d_idx  <- match(d_name, colnames(vcov_mat))
    d_est  <- x$growth_rates[v]
    d_se   <- if (!is.na(d_idx)) sqrt(max(vcov_mat[d_idx, d_idx], 0)) else NA_real_

    rows <- c(rows, list(
      tibble::tibble(
        lineage   = v,
        term      = "intercept",
        estimate  = unname(a_est),
        std.error = a_se,
        conf.low  = if (conf.int) unname(a_est) - z * a_se else NA_real_,
        conf.high = if (conf.int) unname(a_est) + z * a_se else NA_real_
      ),
      tibble::tibble(
        lineage   = v,
        term      = "growth_rate",
        estimate  = unname(d_est),
        std.error = d_se,
        conf.low  = if (conf.int) unname(d_est) - z * d_se else NA_real_,
        conf.high = if (conf.int) unname(d_est) + z * d_se else NA_real_
      )
    ))
  }

  dplyr::bind_rows(rows)
}


#' Glance at an lfq_fit object
#'
#' Returns a single-row tibble of model-level summary statistics.
#'
#' @param x An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return A single-row tibble with columns: `engine`, `n_lineages`,
#'   `n_timepoints`, `nobs`, `df`, `logLik`, `AIC`, `BIC`, `pivot`,
#'   `convergence`.
#'
#' @examples
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
#' fit <- fit_model(sim)
#' glance.lfq_fit(fit)
#'
#' @export
glance.lfq_fit <- function(x, ...) {
  tibble::tibble(
    engine       = x$engine,
    n_lineages   = length(x$lineages),
    n_timepoints = as.integer(x$n_timepoints),
    nobs         = as.integer(x$nobs),
    df           = as.integer(x$df),
    logLik       = x$loglik,
    AIC          = x$aic,
    BIC          = x$bic,
    pivot        = x$pivot,
    convergence  = as.integer(x$convergence)
  )
}


#' Augment data with fitted values from an lfq_fit object
#'
#' @param x An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return A tibble with columns: `.date`, `.lineage`, `.fitted_freq`,
#'   `.observed`, `.pearson_resid`.
#'
#' @examples
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
#' fit <- fit_model(sim)
#' augment.lfq_fit(fit)
#'
#' @export
augment.lfq_fit <- function(x, ...) {
  x$residuals
}
