#' Print an lfq_fit object
#'
#' @param x An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 3, advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
#'   n_timepoints = 15, seed = 42
#' )
#' fit <- fit_model(sim, engine = "mlr")
#' print(fit)
#'
#' @export
print.lfq_fit <- function(x, ...) {
  cli::cli_h2("Lineage frequency model ({x$engine})")

  cli::cli_text(
    "Pivot: {.val {x$pivot}}  |  Lineages: {.val {length(x$lineages)}}  |  Time points: {.val {x$n_timepoints}}"
  )
  cli::cli_text(
    "Date range: {x$date_range[1]} to {x$date_range[2]}"
  )

  cli::cli_h3("Growth rates (per {x$time_scale} days)")

  gr <- sort(x$growth_rates, decreasing = TRUE)
  for (nm in names(gr)) {
    delta <- gr[nm]
    arrow <- if (nm == x$pivot) {
      cli::col_grey("\u25A0")          # filled square for pivot
    } else if (delta > 0) {
      cli::col_green("\u25B2")         # up triangle
    } else {
      cli::col_red("\u25BC")           # down triangle
    }
    label <- if (nm == x$pivot) "(pivot)" else ""
    cli::cli_text(
      "  {arrow} {nm}{label}: {round(delta, 4)}"
    )
  }

  cli::cli_text("")
  cli::cli_text(
    "log-Lik: {round(x$loglik, 1)}  AIC: {round(x$aic, 1)}  BIC: {round(x$bic, 1)}"
  )
  cli::cli_text(
    "Convergence: {if (x$convergence == 0) 'OK' else paste('code', x$convergence)}"
  )

  invisible(x)
}


#' Summarise an lfq_fit object
#'
#' @param object An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `object`.
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 3, advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
#'   n_timepoints = 15, seed = 42
#' )
#' fit <- fit_model(sim, engine = "mlr")
#' summary(fit)
#'
#' @export
summary.lfq_fit <- function(object, ...) {
  print(object)

  cli::cli_h3("Parameter table (non-pivot)")
  non_pivot <- setdiff(object$lineages, object$pivot)
  vcov_m    <- object$vcov_matrix

  for (v in non_pivot) {
    d      <- object$growth_rates[v]
    a      <- object$intercepts[v]
    se_nm  <- paste0("delta_", v)
    se_idx <- match(se_nm, colnames(vcov_m))
    se_d   <- if (!is.na(se_idx)) sqrt(max(vcov_m[se_idx, se_idx], 0)) else NA_real_
    z      <- stats::qnorm((1 + object$ci_level) / 2)
    cli::cli_text(
      "  {v}: delta = {round(d, 4)} (SE {round(se_d, 4)}), alpha = {round(a, 4)}"
    )
  }

  cli::cli_text("")
  cli::cli_text("n obs (sequences): {object$nobs}")

  invisible(object)
}
