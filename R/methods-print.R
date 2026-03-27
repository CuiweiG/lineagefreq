#' Print a lineage frequency model
#'
#' @param x An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
#'   n_timepoints = 15, seed = 42)
#' fit <- fit_model(sim, engine = "mlr")
#' print(fit)
#'
#' @export
print.lfq_fit <- function(x, ...) {
  cat(sprintf("Lineage frequency model (%s)\n", x$engine))
  cli::cli_text(
    "{length(x$lineages)} lineage{?s}, {x$n_timepoints} time point{?s}"
  )
  cli::cli_text("Date range: {x$date_range[1]} to {x$date_range[2]}")
  cli::cli_text("Pivot: {.val {x$pivot}}")
  cat("\n")
  cat(sprintf("Growth rates (per %d-day unit):\n", x$time_scale))

  non_pivot <- setdiff(x$lineages, x$pivot)
  for (v in non_pivot) {
    d     <- x$growth_rates[v]
    arrow <- if (d > 0.01) "\u2191" else if (d < -0.01) "\u2193" else "\u2192"
    cat(sprintf("  %s %s: %s\n", arrow, v, format(d, digits = 4)))
  }

  cat("\n")
  cat(sprintf("AIC: %s; BIC: %s\n",
              format(x$aic, digits = 1),
              format(x$bic, digits = 1)))
  if (x$convergence != 0) {
    cli::cli_alert_warning("Optimizer did not fully converge.")
  }
  invisible(x)
}


#' Summarise a lineage frequency model
#'
#' @param object An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `object`.
#'
#' @examples
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
#'   n_timepoints = 15, seed = 42)
#' fit <- fit_model(sim, engine = "mlr")
#' summary(fit)
#'
#' @export
summary.lfq_fit <- function(object, ...) {
  ga <- growth_advantage(object, type = "growth_rate")
  cat("Lineage Frequency Model Summary\n")
  cat("================================\n")
  cat("Engine:      ", object$engine, "\n")
  cat("Pivot:       ", object$pivot, "\n")
  cat("Lineages:    ", length(object$lineages), "\n")
  cat("Time points: ", object$n_timepoints, "\n")
  cat("Total seqs:  ", object$nobs, "\n")
  cat("Parameters:  ", object$df, "\n")
  cat("Log-lik:     ", format(object$loglik, digits = 4), "\n")
  cat("AIC:         ", format(object$aic, digits = 4), "\n")
  cat("BIC:         ", format(object$bic, digits = 4), "\n")
  cat("\nGrowth rates (per", object$time_scale, "days):\n")
  print(ga, n = Inf)
  invisible(object)
}
