#' Fit a lineage frequency model
#'
#' Unified interface for modeling lineage frequency dynamics. Supports
#' multiple engines that share the same input/output contract.
#'
#' @param data An [lfq_data] object.
#' @param engine Model engine to use:
#'   * `"mlr"` (default): Multinomial logistic regression. Fast,
#'     frequentist, no external dependencies.
#'   * `"hier_mlr"`: Hierarchical MLR with partial pooling across
#'     locations. Requires `.location` column in data.
#'   * `"piantham"`: Piantham approximation converting MLR growth
#'     rates to relative reproduction numbers. Requires
#'     `generation_time` argument.
#'   * `"fga"`: Fixed growth advantage model (Bayesian via CmdStan).
#'     Requires 'CmdStan'; check with [lfq_stan_available()].
#'   * `"garw"`: Growth advantage random walk model (Bayesian via
#'     CmdStan). Allows fitness to change over time.
#' @param pivot Reference lineage name. Growth rates are reported
#'   relative to this lineage (fixed at 0). Default: the lineage with
#'   the highest count at the earliest time point.
#' @param horizon Forecast horizon in days (stored for later use by
#'   [forecast()]). Default 28.
#' @param ci_level Confidence level for intervals. Default 0.95.
#' @param ... Engine-specific arguments passed to the internal engine
#'   function. For `engine = "mlr"`: `window`, `ci_method`,
#'   `laplace_smooth`. For `engine = "piantham"`: `generation_time`
#'   (required). For `engine = "hier_mlr"`: `shrinkage_method`.
#'
#' @return An `lfq_fit` object (S3 class), a list containing:
#' \describe{
#'   \item{engine}{Engine name (character).}
#'   \item{growth_rates}{Named numeric vector of growth rates per
#'     `time_scale` days (pivot = 0).}
#'   \item{intercepts}{Named numeric vector of intercepts.}
#'   \item{pivot}{Name of pivot lineage.}
#'   \item{lineages}{Character vector of all lineage names.}
#'   \item{fitted_values}{Tibble of fitted frequencies.}
#'   \item{residuals}{Tibble with observed, fitted, Pearson residuals.}
#'   \item{vcov_matrix}{Variance-covariance matrix.}
#'   \item{loglik, aic, bic}{Model fit statistics.}
#'   \item{nobs, n_timepoints, df}{Sample and model sizes.}
#'   \item{ci_level, horizon}{As specified.}
#'   \item{call}{The matched call.}
#' }
#'
#' @seealso [growth_advantage()] to extract fitness estimates,
#'   [forecast()] for frequency prediction, [backtest()] for
#'   rolling-origin evaluation.
#'
#' @references
#' Abousamra E, Figgins M, Bedford T (2024). Fitness models provide
#' accurate short-term forecasts of SARS-CoV-2 variant frequency.
#' \emph{PLoS Computational Biology}, 20(9):e1012443.
#' \doi{10.1371/journal.pcbi.1012443}
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 3,
#'   advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
#'   n_timepoints = 15, seed = 42
#' )
#' fit <- fit_model(sim, engine = "mlr")
#' fit
#'
#' @export
fit_model <- function(data,
                      engine   = c("mlr", "hier_mlr", "piantham", "fga", "garw"),
                      pivot    = NULL,
                      horizon  = 28L,
                      ci_level = 0.95,
                      ...) {

  engine <- match.arg(engine)

  if (!is_lfq_data(data)) {
    cli::cli_abort(
      "{.arg data} must be an {.cls lfq_data} object. Use {.fn lfq_data} first."
    )
  }
  assert_prob(ci_level, "ci_level")

  result <- switch(engine,
    mlr      = .engine_mlr(data, pivot = pivot, ci_level = ci_level, ...),
    hier_mlr = .engine_hier_mlr(data, pivot = pivot, ci_level = ci_level, ...),
    piantham = .engine_piantham(data, pivot = pivot, ci_level = ci_level, ...),
    fga      = .engine_fga(data, pivot = pivot, ci_level = ci_level, ...),
    garw     = .engine_garw(data, pivot = pivot, ci_level = ci_level, ...),
    cli::cli_abort("Unknown engine {.val {engine}}.")
  )

  result$engine   <- engine
  result$ci_level <- ci_level
  result$horizon  <- as.integer(horizon)
  result$call     <- match.call()

  structure(result, class = c(paste0("lfq_fit_", engine), "lfq_fit"))
}


#' List available modeling engines
#'
#' Returns information about all modeling engines available in
#' lineagefreq. Core engines (mlr, hier_mlr, piantham) are always
#' available. Bayesian engines (fga, garw) require 'CmdStan'.
#'
#' @return A tibble with columns: `engine`, `type`, `time_varying`,
#'   `available`, `description`.
#'
#' @examples
#' lfq_engines()
#'
#' @export
lfq_engines <- function() {
  stan_ok <- lfq_stan_available()
  tibble::tibble(
    engine       = c("mlr", "hier_mlr", "piantham", "fga", "garw"),
    type         = c("frequentist", "frequentist", "frequentist",
                     "bayesian", "bayesian"),
    time_varying = c(FALSE, FALSE, FALSE, FALSE, TRUE),
    available    = c(TRUE, TRUE, TRUE, stan_ok, stan_ok),
    description  = c(
      "Multinomial logistic regression (MLE)",
      "Hierarchical MLR with empirical Bayes shrinkage",
      "Piantham Rt approximation from MLR growth rates",
      "Fixed growth advantage (Stan)",
      "Growth advantage random walk (Stan)"
    )
  )
}

