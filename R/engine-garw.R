#' GARW engine (internal)
#'
#' Growth Advantage Random Walk model via 'CmdStan'.
#' Allows growth advantages to change over time.
#' Called by fit_model(engine = "garw"). Not exported.
#'
#' @noRd
.engine_garw <- function(data,
                         pivot         = NULL,
                         ci_level      = 0.95,
                         chains        = 4L,
                         iter_warmup   = 1000L,
                         iter_sampling = 1000L,
                         ...) {

  if (!lfq_stan_available()) {
    cli::cli_abort(c(
      "'CmdStan' is required for the GARW engine.",
      "i" = "Install: {.code cmdstanr::install_cmdstan()}"
    ))
  }

  data_fit   <- data[data$.reliable, ]
  stan_input <- .prepare_stan_data(data_fit, pivot)

  mod <- .get_stan_model("garw")

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

  rlang::check_installed("posterior",
    reason = "to extract posterior draws from Stan models."
  )

  lineages   <- stan_input$lineages
  pivot_name <- stan_input$pivot_name
  n_lin      <- stan_input$V
  n_tp       <- stan_input$T
  dates      <- stan_input$dates

  # Time-varying rho: extract final time point as summary
  n_draws <- chains * iter_sampling
  rho_final_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_lin)
  for (v in seq_len(n_lin)) {
    var_name <- paste0("rho[", n_tp, ",", v, "]")
    rho_final_draws[, v] <- as.numeric(
      fit_stan$draws(variables = var_name, format = "matrix")
    )
  }
  colnames(rho_final_draws) <- lineages

  growth_rates <- stats::setNames(
    apply(rho_final_draws, 2, function(x) log(stats::median(x))),
    lineages
  )

  # Time-varying rho trajectory
  alpha_q <- (1 - ci_level) / 2
  rho_trajectory <- list()
  for (t_idx in seq_len(n_tp)) {
    for (v_idx in seq_len(n_lin)) {
      var_name <- paste0("rho[", t_idx, ",", v_idx, "]")
      vals <- as.numeric(
        fit_stan$draws(variables = var_name, format = "matrix")
      )
      rho_trajectory <- c(rho_trajectory, list(tibble::tibble(
        .date      = dates[t_idx],
        .lineage   = lineages[v_idx],
        rho_median = stats::median(vals),
        rho_lower  = stats::quantile(vals, alpha_q,     names = FALSE),
        rho_upper  = stats::quantile(vals, 1 - alpha_q, names = FALSE)
      )))
    }
  }
  rho_trajectory_df <- dplyr::bind_rows(rho_trajectory)

  # sigma_rw posterior
  sigma_draws <- as.numeric(
    fit_stan$draws(variables = "sigma_rw", format = "matrix")
  )
  sigma_rw_summary <- c(
    mean   = mean(sigma_draws),
    median = stats::median(sigma_draws),
    lower  = stats::quantile(sigma_draws, 0.025, names = FALSE),
    upper  = stats::quantile(sigma_draws, 0.975, names = FALSE)
  )

  # Fitted frequencies
  fitted_rows <- list()
  subset_fn <- utils::getFromNamespace("subset_draws", "posterior")
  freq_draws <- fit_stan$draws(variables = "freq_hat",
                               format = "draws_array")
  for (t_idx in seq_len(n_tp)) {
    for (v_idx in seq_len(n_lin)) {
      var_name <- paste0("freq_hat[", t_idx, ",", v_idx, "]")
      vals <- as.numeric(subset_fn(freq_draws, variable = var_name))
      fitted_rows <- c(fitted_rows, list(tibble::tibble(
        .date        = dates[t_idx],
        .lineage     = lineages[v_idx],
        .fitted_freq = stats::median(vals)
      )))
    }
  }
  fitted_df <- dplyr::bind_rows(fitted_rows)

  # Residuals
  obs_df   <- data_fit[, c(".date", ".lineage", ".freq")]
  obs_df   <- dplyr::rename(obs_df, .observed = ".freq")
  resid_df <- dplyr::left_join(fitted_df, obs_df,
                               by = c(".date", ".lineage"))
  resid_df$.pearson_resid <- NA_real_

  # vcov placeholder (diagonal from final rho posterior)
  non_pivot   <- setdiff(lineages, pivot_name)
  n_params    <- 2L * (n_lin - 1L)
  vcov_mat    <- diag(n_params) * 0.01
  param_names <- c(paste0("alpha_", non_pivot),
                   paste0("delta_", non_pivot))

  for (i in seq_along(non_pivot)) {
    v_idx <- match(non_pivot[i], lineages)
    delta_vals <- log(rho_final_draws[, v_idx])
    vcov_mat[length(non_pivot) + i, length(non_pivot) + i] <- stats::var(delta_vals)
  }
  rownames(vcov_mat) <- colnames(vcov_mat) <- param_names

  list(
    growth_rates    = growth_rates,
    intercepts      = stats::setNames(rep(0, n_lin), lineages),
    pivot           = pivot_name,
    lineages        = lineages,
    fitted_values   = fitted_df,
    residuals       = resid_df,
    vcov_matrix     = vcov_mat,
    loglik          = NA_real_,
    aic             = NA_real_,
    bic             = NA_real_,
    nobs            = as.integer(sum(stan_input$Y)),
    n_timepoints    = as.integer(n_tp),
    df              = as.integer(n_params),
    convergence     = 0L,
    ci_method       = "hdi",
    date_range      = range(dates),
    time_scale      = 7,
    stan_fit        = fit_stan,
    rho_trajectory  = rho_trajectory_df,
    sigma_rw        = sigma_rw_summary,
    is_time_varying = TRUE
  )
}
