#' FGA engine (internal)
#'
#' Fixed Growth Advantage model via 'CmdStan'.
#' Called by fit_model(engine = "fga"). Not exported.
#'
#' @noRd
.engine_fga <- function(data,
                        pivot         = NULL,
                        ci_level      = 0.95,
                        chains        = 4L,
                        iter_warmup   = 1000L,
                        iter_sampling = 1000L,
                        ...) {

  if (!lfq_stan_available()) {
    cli::cli_abort(c(
      "'CmdStan' is required for the FGA engine.",
      "i" = "Install: {.code cmdstanr::install_cmdstan()}"
    ))
  }

  # Filter reliable data
  data_fit <- data[data$.reliable, ]

  # Prepare Stan data
  stan_input <- .prepare_stan_data(data_fit, pivot)

  # Compile and sample
  mod <- .get_stan_model("fga")

  fit_stan <- mod$sample(
    data = list(
      T     = stan_input$T,
      V     = stan_input$V,
      Y     = stan_input$Y,
      pivot = stan_input$pivot
    ),
    chains          = chains,
    parallel_chains = min(chains, 4L),
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    refresh         = 0,
    show_messages   = FALSE,
    ...
  )

  # Extract posteriors
  rlang::check_installed("posterior",
    reason = "to extract posterior draws from Stan models."
  )

  lineages   <- stan_input$lineages
  pivot_name <- stan_input$pivot_name
  n_lin      <- stan_input$V

  # rho posteriors
  rho_draws <- fit_stan$draws(variables = "rho", format = "matrix")
  rho_names <- paste0("rho[", seq_len(n_lin), "]")
  colnames(rho_draws) <- rho_names

  rho_summary <- apply(rho_draws, 2, function(x) {
    c(mean   = mean(x),
      median = stats::median(x),
      lower  = stats::quantile(x, (1 - ci_level) / 2, names = FALSE),
      upper  = stats::quantile(x, 1 - (1 - ci_level) / 2, names = FALSE))
  })

  growth_advantages_Rt <- stats::setNames(rho_summary["mean", ], lineages)
  growth_rates <- log(growth_advantages_Rt)

  # Fitted frequencies from posterior
  freq_draws <- fit_stan$draws(variables = "freq_hat",
                               format = "draws_array")

  fitted_rows <- list()
  dates <- stan_input$dates
  for (t_idx in seq_along(dates)) {
    for (v_idx in seq_along(lineages)) {
      var_name <- paste0("freq_hat[", t_idx, ",", v_idx, "]")
      subset_fn <- utils::getFromNamespace("subset_draws", "posterior")
      vals <- as.numeric(subset_fn(freq_draws, variable = var_name))
      fitted_rows <- c(fitted_rows, list(tibble::tibble(
        .date        = dates[t_idx],
        .lineage     = lineages[v_idx],
        .fitted_freq = stats::median(vals)
      )))
    }
  }
  fitted_df <- dplyr::bind_rows(fitted_rows)

  # Build residuals
  obs_df <- data_fit[, c(".date", ".lineage", ".freq")]
  obs_df <- dplyr::rename(obs_df, .observed = ".freq")
  resid_df <- dplyr::left_join(fitted_df, obs_df,
                               by = c(".date", ".lineage"))
  resid_df$.pearson_resid <- NA_real_

  # Build vcov from posterior (diagonal approximation)
  non_pivot  <- setdiff(lineages, pivot_name)
  n_params   <- 2L * (n_lin - 1L)
  vcov_mat   <- diag(n_params) * 0.01
  param_names <- c(paste0("alpha_", non_pivot),
                   paste0("delta_", non_pivot))

  for (i in seq_along(non_pivot)) {
    v     <- non_pivot[i]
    v_idx <- match(v, lineages)
    rho_vals   <- rho_draws[, paste0("rho[", v_idx, "]")]
    delta_vals <- log(rho_vals)
    vcov_mat[length(non_pivot) + i, length(non_pivot) + i] <- stats::var(delta_vals)
  }
  rownames(vcov_mat) <- colnames(vcov_mat) <- param_names

  list(
    growth_rates  = growth_rates,
    intercepts    = stats::setNames(rep(0, n_lin), lineages),
    pivot         = pivot_name,
    lineages      = lineages,
    fitted_values = fitted_df,
    residuals     = resid_df,
    vcov_matrix   = vcov_mat,
    loglik        = NA_real_,
    aic           = NA_real_,
    bic           = NA_real_,
    nobs          = as.integer(sum(stan_input$Y)),
    n_timepoints  = stan_input$n_timepoints,
    df            = as.integer(n_params),
    convergence   = 0L,
    ci_method     = "hdi",
    date_range    = range(dates),
    time_scale    = 7,
    stan_fit      = fit_stan,
    rho_summary   = rho_summary,
    rho_draws     = rho_draws
  )
}
