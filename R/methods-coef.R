#' Extract coefficients from an lfq_fit object
#'
#' Returns all growth rates and intercepts as a named numeric vector.
#'
#' @param object An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return Named numeric vector with elements `alpha_<lineage>` and
#'   `delta_<lineage>` for each non-pivot lineage.
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 3, advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 12, seed = 1
#' )
#' fit <- fit_model(sim, engine = "mlr")
#' coef(fit)
#'
#' @export
coef.lfq_fit <- function(object, ...) {
  non_pivot <- setdiff(object$lineages, object$pivot)
  alphas    <- object$intercepts[non_pivot]
  deltas    <- object$growth_rates[non_pivot]
  names(alphas) <- paste0("alpha_", non_pivot)
  names(deltas) <- paste0("delta_", non_pivot)
  c(alphas, deltas)
}


#' Extract variance-covariance matrix from an lfq_fit object
#'
#' @param object An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return A named numeric matrix.
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 3, advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 12, seed = 1
#' )
#' fit <- fit_model(sim, engine = "mlr")
#' vcov(fit)
#'
#' @export
vcov.lfq_fit <- function(object, ...) {
  object$vcov_matrix
}


#' Extract log-likelihood from an lfq_fit object
#'
#' @param object An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return An object of class `logLik`.
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 3, advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 12, seed = 1
#' )
#' fit <- fit_model(sim, engine = "mlr")
#' logLik(fit)
#'
#' @export
logLik.lfq_fit <- function(object, ...) {
  val <- object$loglik
  attr(val, "df")    <- object$df
  attr(val, "nobs")  <- object$nobs
  class(val) <- "logLik"
  val
}


#' Extract number of observations from an lfq_fit object
#'
#' @param object An [lfq_fit] object.
#' @param ... Ignored.
#'
#' @return Integer scalar (total sequence count used in fitting).
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 3, advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 12, seed = 1
#' )
#' fit <- fit_model(sim, engine = "mlr")
#' nobs(fit)
#'
#' @export
nobs.lfq_fit <- function(object, ...) {
  as.integer(object$nobs)
}
