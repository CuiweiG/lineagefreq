#' Hierarchical MLR engine (internal)
#'
#' Fits multinomial logistic regression independently per location,
#' then applies empirical Bayes shrinkage (DerSimonian-Laird) to
#' pool growth rates across locations. No lme4 dependency.
#'
#' Called by fit_model(engine = "hier_mlr"). Not exported.
#'
#' @noRd
.engine_hier_mlr <- function(data,
                             pivot            = NULL,
                             ci_level         = 0.95,
                             shrinkage_method = "eb",
                             ...) {

  if (!attr(data, "has_location")) {
    cli::cli_abort(
      c("Engine {.val hier_mlr} requires a {.field .location} column.",
        "i" = "Pass {.code location = <col>} when calling {.fn lfq_data}.")
    )
  }

  locations <- sort(unique(data$.location))
  n_loc     <- length(locations)

  if (n_loc < 2L) {
    cli::cli_abort("Hierarchical MLR requires at least 2 locations. Found {n_loc}.")
  }

  lineages <- sort(unique(data$.lineage))

  # --- Step 1: Fit MLR per location ---
  loc_fits <- list()
  for (loc in locations) {
    loc_data <- data[data$.location == loc, ]
    # Re-wrap as lfq_data to preserve class
    loc_data <- structure(
      loc_data,
      class        = class(data),
      lineages     = sort(unique(loc_data$.lineage)),
      date_range   = range(loc_data$.date, na.rm = TRUE),
      n_timepoints = length(unique(loc_data$.date)),
      has_location = TRUE,
      min_total    = attr(data, "min_total")
    )

    fit <- tryCatch(
      .engine_mlr(loc_data, pivot = pivot, ci_level = ci_level, ...),
      error = function(e) {
        cli::cli_warn("MLR failed for location {.val {loc}}: {e$message}")
        NULL
      }
    )
    if (!is.null(fit)) loc_fits[[loc]] <- fit
  }

  if (length(loc_fits) < 2L) {
    cli::cli_abort("Hierarchical MLR needs at least 2 successful location fits.")
  }

  # Use the pivot from the first successful fit
  pivot_used <- loc_fits[[1L]]$pivot
  non_pivot  <- setdiff(lineages, pivot_used)
  ts         <- loc_fits[[1L]]$time_scale

  # --- Step 2: Collect per-location deltas and SEs ---
  # For each non-pivot lineage, gather delta and SE across locations
  loc_names  <- names(loc_fits)

  delta_by_lin <- list()
  se_by_lin    <- list()

  for (v in non_pivot) {
    deltas <- numeric(0)
    ses    <- numeric(0)
    locs   <- character(0)

    for (loc in loc_names) {
      fit_loc <- loc_fits[[loc]]
      if (!v %in% fit_loc$lineages) next

      d      <- fit_loc$growth_rates[v]
      d_name <- paste0("delta_", v)
      d_idx  <- match(d_name, colnames(fit_loc$vcov_matrix))
      se     <- if (!is.na(d_idx)) sqrt(max(fit_loc$vcov_matrix[d_idx, d_idx], 0)) else NA_real_

      if (!is.na(se) && se > 0) {
        deltas <- c(deltas, d)
        ses    <- c(ses, se)
        locs   <- c(locs, loc)
      }
    }
    delta_by_lin[[v]] <- list(delta = deltas, se = ses, loc = locs)
  }

  # --- Step 3: DerSimonian-Laird meta-analysis per lineage ---
  global_delta    <- stats::setNames(rep(0, length(lineages)), lineages)
  global_se       <- stats::setNames(rep(NA_real_, length(lineages)), lineages)
  shrunk_by_loc   <- list()

  for (v in non_pivot) {
    info <- delta_by_lin[[v]]
    if (length(info$delta) < 2L) {
      # Not enough locations; use whatever we have
      if (length(info$delta) == 1L) {
        global_delta[v] <- info$delta[1]
        global_se[v]    <- info$se[1]
      }
      next
    }

    d_vec  <- info$delta
    se_vec <- info$se
    w_vec  <- 1 / se_vec^2

    # Weighted mean (fixed-effect)
    d_fe <- sum(w_vec * d_vec) / sum(w_vec)

    # DerSimonian-Laird tau^2
    Q      <- sum(w_vec * (d_vec - d_fe)^2)
    k      <- length(d_vec)
    c_dl   <- sum(w_vec) - sum(w_vec^2) / sum(w_vec)
    tau_sq <- max((Q - (k - 1)) / c_dl, 0)

    # Random-effects weights
    w_re    <- 1 / (se_vec^2 + tau_sq)
    d_re    <- sum(w_re * d_vec) / sum(w_re)
    se_re   <- sqrt(1 / sum(w_re))

    global_delta[v] <- d_re
    global_se[v]    <- se_re

    # Shrinkage per location
    for (j in seq_along(info$loc)) {
      loc <- info$loc[j]
      lambda_j <- se_vec[j]^2 / (se_vec[j]^2 + tau_sq)
      d_shrunk <- (1 - lambda_j) * d_vec[j] + lambda_j * d_re

      if (is.null(shrunk_by_loc[[loc]])) {
        shrunk_by_loc[[loc]] <- stats::setNames(rep(0, length(lineages)), lineages)
      }
      shrunk_by_loc[[loc]][v] <- d_shrunk
    }
  }

  # --- Step 4: Build vcov for global estimates ---
  n_params    <- 2L * (length(non_pivot))
  param_names <- c(paste0("alpha_", non_pivot), paste0("delta_", non_pivot))
  vcov_mat    <- matrix(0, nrow = n_params, ncol = n_params,
                        dimnames = list(param_names, param_names))

  for (i in seq_along(non_pivot)) {
    v      <- non_pivot[i]
    se_val <- global_se[v]
    if (!is.na(se_val)) {
      d_idx <- length(non_pivot) + i
      vcov_mat[d_idx, d_idx] <- se_val^2
    }
  }

  # --- Step 5: Build intercepts (from global weighted average) ---
  # Collect intercepts similarly
  global_alpha <- stats::setNames(rep(0, length(lineages)), lineages)
  for (v in non_pivot) {
    alphas <- numeric(0)
    wts    <- numeric(0)
    for (loc in loc_names) {
      fit_loc <- loc_fits[[loc]]
      if (v %in% fit_loc$lineages) {
        a      <- fit_loc$intercepts[v]
        a_name <- paste0("alpha_", v)
        a_idx  <- match(a_name, colnames(fit_loc$vcov_matrix))
        se     <- if (!is.na(a_idx)) sqrt(max(fit_loc$vcov_matrix[a_idx, a_idx], 0)) else NA_real_
        if (!is.na(se) && se > 0) {
          alphas <- c(alphas, a)
          wts    <- c(wts, 1 / se^2)
        }
      }
    }
    if (length(alphas) > 0) {
      global_alpha[v] <- sum(wts * alphas) / sum(wts)
      a_se <- sqrt(1 / sum(wts))
      a_idx_vc <- match(paste0("alpha_", v), param_names)
      vcov_mat[a_idx_vc, a_idx_vc] <- a_se^2
    }
  }

  # --- Step 6: Fitted values (using global params) ---
  all_dates <- sort(unique(data$.date))
  t_vec     <- as.numeric(all_dates - min(all_dates)) / ts

  fitted_rows <- list()
  for (i in seq_along(all_dates)) {
    log_num   <- global_alpha + global_delta * t_vec[i]
    log_denom <- log_sum_exp(log_num)
    probs     <- exp(log_num - log_denom)
    for (v in lineages) {
      fitted_rows <- c(fitted_rows, list(tibble::tibble(
        .date        = all_dates[i],
        .lineage     = v,
        .fitted_freq = unname(probs[v])
      )))
    }
  }
  fitted_df <- dplyr::bind_rows(fitted_rows)

  # --- Step 7: Residuals ---
  obs_df <- data[, c(".date", ".lineage", ".freq")]
  obs_df <- dplyr::rename(obs_df, .observed = ".freq")
  # Average observed freq across locations for global residuals
  obs_agg <- obs_df |>
    dplyr::group_by(.data$.date, .data$.lineage) |>
    dplyr::summarise(.observed = mean(.data$.observed, na.rm = TRUE),
                     .groups = "drop")

  resid_df <- dplyr::left_join(fitted_df, obs_agg, by = c(".date", ".lineage"))
  resid_df$.pearson_resid <- NA_real_

  # --- Step 8: Model statistics ---
  # Aggregate loglik from per-location fits
  total_loglik <- sum(vapply(loc_fits, function(f) f$loglik, numeric(1)))
  total_nobs   <- sum(vapply(loc_fits, function(f) f$nobs, integer(1)))
  n_tp         <- length(all_dates)
  aic_val      <- 2 * n_params - 2 * total_loglik
  bic_val      <- log(n_tp * n_loc) * n_params - 2 * total_loglik

  # --- Return ---
  list(
    growth_rates    = global_delta,
    intercepts      = global_alpha,
    pivot           = pivot_used,
    lineages        = lineages,
    fitted_values   = fitted_df,
    residuals       = resid_df,
    vcov_matrix     = vcov_mat,
    loglik          = total_loglik,
    aic             = aic_val,
    bic             = bic_val,
    nobs            = as.integer(total_nobs),
    n_timepoints    = as.integer(n_tp),
    df              = as.integer(n_params),
    convergence     = 0L,
    ci_method       = "wald",
    date_range      = range(all_dates),
    time_scale      = ts,
    locations        = locations,
    n_locations      = n_loc,
    location_fits    = loc_fits,
    shrunk_by_loc    = shrunk_by_loc,
    shrinkage_method = shrinkage_method
  )
}
