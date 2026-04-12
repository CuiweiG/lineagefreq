#' Fit model with Deep Mutational Scanning priors
#'
#' A penalised multinomial logistic regression engine that
#' incorporates Deep Mutational Scanning (DMS) escape scores as
#' informative priors on variant fitness. This is valuable for
#' early-emergence scenarios where a new lineage has few observed
#' sequences but laboratory-measured phenotypic data (e.g., ACE2
#' binding affinity, antibody escape) are available.
#'
#' The approach uses penalised maximum likelihood where the penalty
#' is proportional to the squared difference between the estimated
#' growth rate and the DMS-derived prior. This implements an
#' empirical Bayes shrinkage: with abundant data, the penalty has
#' little effect; with sparse data, estimates are pulled toward the
#' DMS prior.
#'
#' @param data An \code{lfq_data} object.
#' @param dms_scores Named numeric vector of DMS-derived fitness
#'   priors. Names correspond to lineage identifiers. Values are
#'   on the log growth rate scale (positive = fitter than average).
#'   Lineages not in the vector receive a prior of 0.
#' @param lambda Regularisation strength (penalty weight). Default
#'   1.0. Larger values pull estimates more strongly toward the DMS
#'   prior. At \code{lambda = 0}, the result is identical to the
#'   standard MLR engine.
#' @param pivot Reference lineage name. Default \code{NULL}
#'   (automatic selection).
#' @param ci_level Confidence level. Default 0.95.
#'
#' @return An \code{lfq_fit} object compatible with all downstream
#'   functions (\code{forecast}, \code{growth_advantage}, etc.).
#'
#' @details
#' The penalised log-likelihood is:
#' \deqn{\ell_{\text{pen}}(\alpha, \delta) = \ell(\alpha, \delta) - \frac{\lambda}{2} \sum_v (\delta_v - \mu_v)^2}
#' where \eqn{\ell} is the standard multinomial log-likelihood,
#' \eqn{\delta_v} is the growth rate for lineage \eqn{v}, and
#' \eqn{\mu_v} is the DMS prior. The Hessian is adjusted
#' accordingly, ensuring correct confidence interval widths.
#'
#' @references
#' Dadonaite B, Crawford KHD, Radford CE, et al. (2023). A
#' pseudovirus system enables deep mutational scanning of the full
#' SARS-CoV-2 spike. \emph{Cell}, 186(6), 1263--1278.
#' \doi{10.1016/j.cell.2023.02.001}
#'
#' Bloom JD, Neher RA (2023). Fitness effects of mutations to
#' SARS-CoV-2 proteins. \emph{Virus Evolution}, 9(2), vead055.
#' \doi{10.1093/ve/vead055}
#'
#' @seealso \code{\link{fit_model}} for the standard MLR engine.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.3, "B" = 0.9),
#'   n_timepoints = 8, total_per_tp = 100, seed = 1)
#' # DMS suggests lineage A has fitness advantage
#' dms <- c("A" = 0.04, "B" = -0.02)
#' fit_dms <- fit_dms_prior(sim, dms_scores = dms, lambda = 2)
#' growth_advantage(fit_dms)
#' }
#'
#' @export
fit_dms_prior <- function(data,
                          dms_scores,
                          lambda   = 1.0,
                          pivot    = NULL,
                          ci_level = 0.95) {

  if (!is_lfq_data(data)) {
    cli::cli_abort("{.arg data} must be an {.cls lfq_data} object.")
  }
  if (!is.numeric(dms_scores) || is.null(names(dms_scores))) {
    cli::cli_abort("{.arg dms_scores} must be a named numeric vector.")
  }
  assert_prob(ci_level, "ci_level")

  lineages <- attr(data, "lineages")

  # Select pivot
  if (is.null(pivot)) {
    first_date <- min(data$.date)
    first_obs  <- data[data$.date == first_date, ]
    pivot <- first_obs$.lineage[which.max(first_obs$.count)]
  }

  non_pivot <- setdiff(lineages, pivot)
  n_lin     <- length(lineages)
  n_par     <- 2L * (n_lin - 1L)

  # Build count matrix and time vector
  dates_sorted <- sort(unique(data$.date))
  t0 <- min(dates_sorted)
  ts <- 7L

  count_mat <- matrix(0L, nrow = length(dates_sorted), ncol = n_lin,
                      dimnames = list(NULL, lineages))
  for (i in seq_along(dates_sorted)) {
    d <- dates_sorted[i]
    sub <- data[data$.date == d, ]
    for (lin in lineages) {
      row <- sub[sub$.lineage == lin, ]
      if (nrow(row) == 1L) count_mat[i, lin] <- row$.count
    }
  }

  t_vec <- as.numeric(dates_sorted - t0) / ts
  row_totals <- rowSums(count_mat)

  # DMS priors for non-pivot lineages
  dms_prior <- vapply(non_pivot, function(lin) {
    if (lin %in% names(dms_scores)) dms_scores[lin] else 0
  }, numeric(1L))

  # Penalised negative log-likelihood
  nll_pen <- function(par) {
    alphas <- c(stats::setNames(0, pivot),
                stats::setNames(par[seq_len(n_lin - 1L)], non_pivot))
    deltas <- c(stats::setNames(0, pivot),
                stats::setNames(par[(n_lin - 1L) + seq_len(n_lin - 1L)],
                                non_pivot))

    nll <- 0
    for (i in seq_along(t_vec)) {
      log_num <- alphas[lineages] + deltas[lineages] * t_vec[i]
      log_denom <- log_sum_exp(log_num)
      log_probs <- log_num - log_denom
      nll <- nll - sum(count_mat[i, ] * log_probs[lineages])
    }

    # Penalty: lambda/2 * sum((delta_v - mu_v)^2)
    delta_np <- par[(n_lin - 1L) + seq_len(n_lin - 1L)]
    penalty <- (lambda / 2) * sum((delta_np - dms_prior)^2)

    nll + penalty
  }

  # Initial parameters
  init <- rep(0, n_par)

  opt <- stats::optim(init, nll_pen, method = "BFGS", hessian = TRUE)

  # Extract parameters
  alpha_est <- c(stats::setNames(0, pivot),
                 stats::setNames(opt$par[seq_len(n_lin - 1L)], non_pivot))
  delta_est <- c(stats::setNames(0, pivot),
                 stats::setNames(opt$par[(n_lin - 1L) + seq_len(n_lin - 1L)],
                                 non_pivot))

  # Variance-covariance from penalised Hessian
  vcov_mat <- tryCatch({
    solve(opt$hessian)
  }, error = function(e) {
    MASS::ginv(opt$hessian)
  })

  param_names <- c(paste0("alpha_", non_pivot), paste0("delta_", non_pivot))
  dimnames(vcov_mat) <- list(param_names, param_names)

  # Fitted values
  fitted_rows <- list()
  resid_rows  <- list()
  for (i in seq_along(dates_sorted)) {
    log_num <- alpha_est[lineages] + delta_est[lineages] * t_vec[i]
    log_denom <- log_sum_exp(log_num)
    fitted_freq <- exp(log_num - log_denom)

    for (lin in lineages) {
      obs_freq <- count_mat[i, lin] / max(row_totals[i], 1L)
      ff <- unname(fitted_freq[lin])
      fitted_rows <- c(fitted_rows, list(tibble::tibble(
        .date = dates_sorted[i], .lineage = lin, .fitted_freq = ff
      )))
      pres <- if (ff > 0 && ff < 1) {
        (obs_freq - ff) / sqrt(ff * (1 - ff) / max(row_totals[i], 1L))
      } else 0
      resid_rows <- c(resid_rows, list(tibble::tibble(
        .date = dates_sorted[i], .lineage = lin,
        .observed = obs_freq, .fitted_freq = ff, .pearson_resid = pres
      )))
    }
  }

  loglik <- -opt$value + (lambda / 2) * sum(
    (opt$par[(n_lin - 1L) + seq_len(n_lin - 1L)] - dms_prior)^2)

  result <- list(
    engine        = "dms_prior",
    growth_rates  = delta_est[lineages],
    intercepts    = alpha_est[lineages],
    pivot         = pivot,
    lineages      = lineages,
    fitted_values = dplyr::bind_rows(fitted_rows),
    residuals     = dplyr::bind_rows(resid_rows),
    vcov_matrix   = vcov_mat,
    loglik        = loglik,
    aic           = 2 * n_par - 2 * loglik,
    bic           = log(sum(row_totals)) * n_par - 2 * loglik,
    nobs          = as.integer(sum(row_totals)),
    n_timepoints  = length(dates_sorted),
    df            = as.integer(n_par),
    convergence   = opt$convergence,
    ci_method     = "wald",
    date_range    = range(dates_sorted),
    time_scale    = ts,
    ci_level      = ci_level,
    horizon       = 28L,
    call          = match.call(),
    dms_scores    = dms_scores,
    lambda        = lambda
  )

  structure(result, class = c("lfq_fit_dms_prior", "lfq_fit"))
}
