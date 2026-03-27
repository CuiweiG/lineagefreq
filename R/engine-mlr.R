#' MLR engine (internal)
#'
#' Fits multinomial logistic regression to lineage frequency data.
#' Called by fit_model(engine = "mlr"). Not exported.
#'
#' @noRd
.engine_mlr <- function(data,
                        pivot      = NULL,
                        ci_level   = 0.95,
                        window     = NULL,
                        ci_method  = "wald",
                        laplace_smooth = 0) {

  # --- Apply time window ---
  if (!is.null(window)) {
    assert_pos_int(window, "window")
    max_date <- max(data$.date)
    data <- data[data$.date >= max_date - window, ]
  }

  # --- Filter to reliable time points ---
  data_fit <- data[data$.reliable, ]
  if (nrow(data_fit) == 0) {
    cli::cli_abort("No reliable time points after filtering.")
  }

  lineages <- sort(unique(data_fit$.lineage))
  n_lin    <- length(lineages)
  dates    <- sort(unique(data_fit$.date))
  n_tp     <- length(dates)

  if (n_lin < 2) cli::cli_abort("At least 2 lineages required.")
  if (n_tp  < 2) cli::cli_abort("At least 2 time points required.")

  # --- Select pivot ---
  # Use only the FIRST time point to avoid bias from lineages that
  # have already grown substantially by the end of the first quarter.
  if (is.null(pivot)) {
    first_date    <- dates[1L]
    first_tp_data <- data_fit[data_fit$.date == first_date, ]
    totals_by_lin <- tapply(first_tp_data$.count, first_tp_data$.lineage,
                            sum, na.rm = TRUE)
    pivot <- names(which.max(totals_by_lin))
  }
  if (!pivot %in% lineages) {
    cli::cli_abort(
      "Pivot {.val {pivot}} not found. Available: {.val {lineages}}."
    )
  }

  non_pivot <- setdiff(lineages, pivot)
  n_params  <- 2L * (n_lin - 1L)

  # --- Reshape to wide count matrix ---
  wide <- data_fit[, c(".date", ".lineage", ".count")]
  wide <- tidyr::pivot_wider(
    wide,
    names_from  = ".lineage",
    values_from = ".count",
    values_fill = 0L
  )
  wide <- dplyr::arrange(wide, .data$.date)

  count_mat  <- as.matrix(wide[, lineages, drop = FALSE])
  row_totals <- rowSums(count_mat)

  # Apply Laplace smoothing if requested
  if (laplace_smooth > 0) {
    count_mat  <- count_mat + laplace_smooth
    row_totals <- rowSums(count_mat)
  }

  # --- Time variable (in weeks) ---
  time_scale <- 7
  t_vec <- as.numeric(wide$.date - min(wide$.date)) / time_scale

  # --- Negative log-likelihood ---
  nll <- function(par) {
    alphas <- par[seq_len(n_lin - 1L)]
    deltas <- par[(n_lin - 1L) + seq_len(n_lin - 1L)]

    alpha_full <- stats::setNames(c(0, alphas), c(pivot, non_pivot))[lineages]
    delta_full <- stats::setNames(c(0, deltas), c(pivot, non_pivot))[lineages]

    ll <- 0
    for (i in seq_len(n_tp)) {
      log_num   <- alpha_full + delta_full * t_vec[i]
      log_denom <- log_sum_exp(log_num)
      log_probs <- log_num - log_denom
      ll <- ll + sum(count_mat[i, ] * log_probs)
    }
    -ll
  }

  # --- Initial values ---
  mean_freq  <- colMeans(count_mat / row_totals)
  mean_freq  <- pmax(mean_freq, 1e-6)
  init_alpha <- log(mean_freq[non_pivot] / mean_freq[pivot])
  init_delta <- rep(0, n_lin - 1L)
  init_par   <- c(unname(init_alpha), init_delta)

  # --- Optimize ---
  opt <- stats::optim(
    par     = init_par,
    fn      = nll,
    method  = "BFGS",
    hessian = TRUE,
    control = list(maxit = 5000, reltol = 1e-12)
  )

  if (opt$convergence != 0) {
    cli::cli_warn("Optimizer did not converge (code {opt$convergence}).")
  }

  # --- Extract estimates ---
  alpha_est <- opt$par[seq_len(n_lin - 1L)]
  delta_est <- opt$par[(n_lin - 1L) + seq_len(n_lin - 1L)]
  names(alpha_est) <- non_pivot
  names(delta_est) <- non_pivot

  growth_rates <- stats::setNames(
    c(0, delta_est), c(pivot, non_pivot)
  )[lineages]

  intercepts <- stats::setNames(
    c(0, alpha_est), c(pivot, non_pivot)
  )[lineages]

  # --- Variance-covariance matrix ---
  hess <- opt$hessian
  if (any(is.na(hess)) || rcond(hess) < .Machine$double.eps) {
    hess <- numDeriv::hessian(nll, opt$par)
  }
  vcov_mat <- tryCatch(
    solve(hess),
    error = function(e) {
      cli::cli_warn("Hessian singular; using generalized inverse.")
      MASS::ginv(hess)
    }
  )
  param_names <- c(paste0("alpha_", non_pivot),
                   paste0("delta_", non_pivot))
  rownames(vcov_mat) <- colnames(vcov_mat) <- param_names

  # --- Fitted values ---
  fitted_rows <- list()
  for (i in seq_len(n_tp)) {
    log_num   <- intercepts + growth_rates * t_vec[i]
    log_denom <- log_sum_exp(log_num)
    probs     <- exp(log_num - log_denom)
    for (v in lineages) {
      fitted_rows <- c(fitted_rows, list(tibble::tibble(
        .date        = dates[i],
        .lineage     = v,
        .fitted_freq = unname(probs[v])
      )))
    }
  }
  fitted_df <- dplyr::bind_rows(fitted_rows)

  # --- Residuals ---
  obs_df <- data_fit[, c(".date", ".lineage", ".freq")]
  obs_df <- dplyr::rename(obs_df, .observed = ".freq")

  resid_df <- dplyr::left_join(fitted_df, obs_df,
                               by = c(".date", ".lineage"))

  # Pearson residuals
  resid_df$.pearson_resid <- NA_real_
  for (i in seq_len(nrow(resid_df))) {
    d     <- resid_df$.date[i]
    t_idx <- match(d, dates)
    if (!is.na(t_idx) && !is.na(resid_df$.observed[i])) {
      p_hat <- resid_df$.fitted_freq[i]
      n_t   <- row_totals[t_idx]
      denom <- sqrt(p_hat * (1 - p_hat) / n_t)
      if (denom > 0) {
        resid_df$.pearson_resid[i] <-
          (resid_df$.observed[i] - p_hat) / denom
      }
    }
  }

  # --- Model fit statistics ---
  loglik_val <- -opt$value
  aic_val    <- 2 * n_params - 2 * loglik_val
  bic_val    <- log(n_tp) * n_params - 2 * loglik_val

  # --- Return list ---
  list(
    growth_rates  = growth_rates,
    intercepts    = intercepts,
    pivot         = pivot,
    lineages      = lineages,
    fitted_values = fitted_df,
    residuals     = resid_df,
    vcov_matrix   = vcov_mat,
    loglik        = loglik_val,
    aic           = aic_val,
    bic           = bic_val,
    nobs          = as.integer(sum(count_mat)),
    n_timepoints  = as.integer(n_tp),
    df            = as.integer(n_params),
    convergence   = opt$convergence,
    ci_method     = ci_method,
    date_range    = range(dates),
    time_scale    = time_scale,
    .optim        = opt,
    .count_mat    = count_mat,
    .t_vec        = t_vec,
    .row_totals   = row_totals
  )
}
