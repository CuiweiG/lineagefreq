#' Extract coefficients from a lineage frequency model
#'
#' @param object An `lfq_fit` object.
#' @param type What to return: `"growth_rate"` (default) for growth
#'   rates only, or `"all"` for intercepts and growth rates.
#' @param ... Ignored.
#'
#' @return Named numeric vector.
#'
#' @examples
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
#' fit <- fit_model(sim)
#' coef(fit)
#'
#' @export
coef.lfq_fit <- function(object, type = c("growth_rate", "all"), ...) {
  type <- match.arg(type)
  if (type == "growth_rate") {
    return(object$growth_rates)
  }
  c(object$intercepts, object$growth_rates)
}


#' @return A named numeric matrix.
#' @export
vcov.lfq_fit <- function(object, ...) {
  object$vcov_matrix
}


#' @return An object of class `logLik`.
#' @export
logLik.lfq_fit <- function(object, ...) {
  ll <- object$loglik
  attr(ll, "df")   <- object$df
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  ll
}


#' @return Integer scalar (total sequence count).
#' @export
nobs.lfq_fit <- function(object, ...) {
  object$nobs
}
