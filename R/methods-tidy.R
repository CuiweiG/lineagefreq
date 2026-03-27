#' Tidy an lfq_fit object into a parameter table
#'
#' Returns a tidy tibble of parameter estimates with standard errors
#' and confidence intervals.
#'
#' @param x An [lfq_fit] object.
#' @param conf.int Logical; include confidence intervals? Default `TRUE`.
#' @param conf.level Confidence level. Default uses `x$ci_level`.
#' @param ... Ignored.
#'
#' @return A tibble with columns: `lineage`, `term`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, and (if `conf.int = TRUE`)
#'   `conf.low`, `conf.high`.
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 3, advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 12, seed = 1
#' )
#' fit <- fit_model(sim, engine = "mlr")
#' tidy.lfq_fit(fit)
#'
#' @export
tidy.lfq_fit <- function(x, conf.int = TRUE,
                         conf.level = NULL, ...) {
  conf.level <- conf.level %||% x$ci_level
  z          <- stats::qnorm((1 + conf.level) / 2)
  non_pivot  <- setdiff(x$lineages, x$pivot)
  vcov_m     <- x$vcov_matrix

  rows <- list()
  for (v in non_pivot) {
    for (term_type in c("alpha", "delta")) {
      nm    <- paste0(term_type, "_", v)
      est   <- if (term_type == "alpha") x$intercepts[v] else x$growth_rates[v]
      idx   <- match(nm, colnames(vcov_m))
      se    <- if (!is.na(idx)) sqrt(max(vcov_m[idx, idx], 0)) else NA_real_
      stat  <- if (!is.na(se) && se > 0) est / se else NA_real_
      pval  <- if (!is.na(stat)) 2 * stats::pnorm(-abs(stat)) else NA_real_
      row   <- tibble::tibble(
        lineage   = v,
        term      = nm,
        estimate  = unname(est),
        std.error = se,
        statistic = stat,
        p.value   = pval
      )
      if (conf.int) {
        row$conf.low  <- unname(est) - z * se
        row$conf.high <- unname(est) + z * se
      }
      rows <- c(rows, list(row))
    }
  }
  dplyr::bind_rows(rows)
}


#' Glance at an lfq_fit object
#'
#' Returns a one-row summary tibble of model-level statistics.
#'
#' @param x An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return A one-row tibble with columns: `engine`, `n_lineages`,
#'   `n_timepoints`, `nobs`, `df`, `logLik`, `AIC`, `BIC`,
#'   `convergence`.
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 3, advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 12, seed = 1
#' )
#' fit <- fit_model(sim, engine = "mlr")
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
    convergence  = as.integer(x$convergence)
  )
}


#' Augment an lfq_fit object with residuals
#'
#' Returns the residuals tibble from a fitted model, with fitted
#' values and Pearson residuals attached.
#'
#' @param x An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return A tibble with columns `.date`, `.lineage`, `.fitted_freq`,
#'   `.observed`, `.pearson_resid`.
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 3, advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 12, seed = 1
#' )
#' fit <- fit_model(sim, engine = "mlr")
#' augment.lfq_fit(fit)
#'
#' @export
augment.lfq_fit <- function(x, ...) {
  x$residuals
}
