#' Extract growth advantage estimates
#'
#' Computes relative fitness of each lineage from a fitted model.
#' Supports four output formats for different use cases.
#'
#' @param fit An `lfq_fit` object returned by [fit_model()].
#' @param type Output type:
#'   * `"growth_rate"` (default): raw growth rate delta per
#'     `time_scale` days (typically per week).
#'   * `"relative_Rt"`: relative effective reproduction number.
#'     Requires `generation_time`.
#'   * `"selection_coefficient"`: relative Rt minus 1.
#'     Requires `generation_time`.
#'   * `"doubling_time"`: days for frequency ratio vs pivot to
#'     double (positive = growing) or halve (negative = declining).
#' @param generation_time Mean generation time in days. Required for
#'   `type = "relative_Rt"` and `"selection_coefficient"`. Common
#'   values: SARS-CoV-2 approximately 5 days, influenza approximately
#'   3 days.
#' @param ci_level Confidence level for intervals. Default uses the
#'   level from the fitted model.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{lineage}{Lineage name.}
#'   \item{estimate}{Point estimate.}
#'   \item{lower}{Lower confidence bound.}
#'   \item{upper}{Upper confidence bound.}
#'   \item{type}{Type of estimate.}
#'   \item{pivot}{Name of pivot (reference) lineage.}
#' }
#'
#' @details
#' The `"relative_Rt"` and `"selection_coefficient"` types use the
#' Piantham approximation (Piantham et al. 2022,
#' \doi{10.3390/v14112556}), which assumes:
#' \enumerate{
#'   \item Variants are in their exponential growth phase (not
#'     saturating due to population-level immunity).
#'   \item All variants share the same generation time distribution.
#'   \item Growth advantage reflects transmissibility differences;
#'     immune escape is not separately identified.
#' }
#' Confidence intervals for `"relative_Rt"` are computed in
#' log-space (Wald intervals on log-Rt), which is more accurate
#' than the linear delta method for ratios.
#'
#' @references
#' Piantham C, Linton NM, Nishiura H (2022). Predicting the
#' trajectory of replacements of SARS-CoV-2 variants using
#' relative reproduction numbers. \emph{Viruses}, 14(11):2556.
#' \doi{10.3390/v14112556}
#'
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
#'
#' # Growth rates per week
#' growth_advantage(fit, type = "growth_rate")
#'
#' # Relative Rt (needs generation time)
#' growth_advantage(fit, type = "relative_Rt", generation_time = 5)
#'
#' @export
growth_advantage <- function(fit,
                             type = c("growth_rate", "relative_Rt",
                                      "selection_coefficient",
                                      "doubling_time"),
                             generation_time = NULL,
                             ci_level = NULL) {

  if (!inherits(fit, "lfq_fit")) {
    cli::cli_abort("{.arg fit} must be an {.cls lfq_fit} object.")
  }

  type     <- match.arg(type)
  ci_level <- ci_level %||% fit$ci_level

  if (type %in% c("relative_Rt", "selection_coefficient") &&
      is.null(generation_time)) {
    cli::cli_abort(
      "{.arg generation_time} is required for type {.val {type}}."
    )
  }

  delta     <- fit$growth_rates
  pivot     <- fit$pivot
  ts        <- fit$time_scale
  vcov_m    <- fit$vcov_matrix
  non_pivot <- setdiff(fit$lineages, pivot)

  # Standard errors for delta (non-pivot only)
  delta_names <- paste0("delta_", non_pivot)
  delta_se <- stats::setNames(
    vapply(delta_names, function(nm) {
      idx <- match(nm, colnames(vcov_m))
      if (!is.na(idx)) sqrt(max(vcov_m[idx, idx], 0)) else NA_real_
    }, numeric(1L)),
    non_pivot
  )

  z <- stats::qnorm((1 + ci_level) / 2)

  # Build result row by row
  rows <- vector("list", length(fit$lineages))

  for (i in seq_along(fit$lineages)) {
    v <- fit$lineages[i]
    d <- delta[v]

    if (v == pivot) {
      ref_val <- switch(type,
        growth_rate           = 0,
        relative_Rt           = 1,
        selection_coefficient = 0,
        doubling_time         = Inf
      )
      rows[[i]] <- tibble::tibble(
        lineage  = v,
        estimate = ref_val,
        lower    = ref_val,
        upper    = ref_val,
        type     = type,
        pivot    = pivot
      )
      next
    }

    se <- delta_se[v]

    rows[[i]] <- switch(type,

      growth_rate = tibble::tibble(
        lineage  = v,
        estimate = d,
        lower    = d - z * se,
        upper    = d + z * se,
        type     = type,
        pivot    = pivot
      ),

      relative_Rt = {
        tau <- generation_time / ts
        rho <- exp(d * tau)
        tibble::tibble(
          lineage  = v,
          estimate = rho,
          lower    = exp((d - z * se) * tau),
          upper    = exp((d + z * se) * tau),
          type     = type,
          pivot    = pivot
        )
      },

      selection_coefficient = {
        tau <- generation_time / ts
        tibble::tibble(
          lineage  = v,
          estimate = exp(d * tau) - 1,
          lower    = exp((d - z * se) * tau) - 1,
          upper    = exp((d + z * se) * tau) - 1,
          type     = type,
          pivot    = pivot
        )
      },

      doubling_time = {
        if (abs(d) < 1e-10) {
          tibble::tibble(
            lineage  = v,
            estimate = Inf,
            lower    = Inf,
            upper    = Inf,
            type     = type,
            pivot    = pivot
          )
        } else {
          dt_est    <- log(2) / abs(d) * ts * sign(d)
          d_lo      <- d - z * se
          d_hi      <- d + z * se
          dt_lo <- if (abs(d_lo) > 1e-10) log(2) / abs(d_lo) * ts * sign(d_lo) else Inf
          dt_hi <- if (abs(d_hi) > 1e-10) log(2) / abs(d_hi) * ts * sign(d_hi) else Inf
          bounds <- sort(c(dt_lo, dt_hi))
          tibble::tibble(
            lineage  = v,
            estimate = dt_est,
            lower    = bounds[1],
            upper    = bounds[2],
            type     = type,
            pivot    = pivot
          )
        }
      }
    )
  }

  dplyr::bind_rows(rows)
}
