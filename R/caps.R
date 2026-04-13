#' Compositional Adaptive Prediction Sets (CAPS)
#'
#' Produces calibrated prediction intervals for variant frequencies using
#' a unified framework that combines horizon-specific variance correction,
#' temporal adaptation, and compositional dimension adjustment.
#'
#' CAPS addresses three limitations of existing conformal methods:
#' (1) standard conformal ignores horizon-dependent miscalibration;
#' (2) ACI uses a fixed learning rate across all horizons;
#' (3) marginal conformal violates the compositional constraint.
#'
#' The prediction radius is:
#'   \eqn{r(t, h) = r_{base}(h) \times \phi(t, h) \times \psi(K)}
#' where \eqn{r_{base}} is initialised from the variance ratio \eqn{R(h)},
#' \eqn{\phi} is updated via horizon-specific ACI, and \eqn{\psi} corrects
#' for the dimensionality of the ILR projection.
#'
#' @param fit An \code{lfq_fit} object from \code{fit_model()}.
#' @param data An \code{lfq_data} object (the data used to fit the model).
#' @param horizons Integer vector of forecast horizons in days.
#' @param alpha Nominal miscoverage rate (default 0.05 for 95\% coverage).
#' @param gamma_0 Base ACI learning rate (default 0.05).
#' @param method One of \code{"caps"} (full framework),
#'   \code{"caps_static"} (no temporal adaptation), or
#'   \code{"caps_marginal"} (no compositional correction).
#' @param n_sim Number of parametric simulations for base intervals
#'   (default 1000).
#' @param ... Passed to \code{backtest()}.
#'
#' @return An object of class \code{caps_forecast} containing:
#'   \describe{
#'     \item{forecasts}{Tibble of point forecasts with CAPS intervals.}
#'     \item{params}{Per-horizon parameters (R_hat, r_base, gamma_h).}
#'     \item{R_hat}{Variance ratios by horizon.}
#'     \item{psi_K}{Compositional dimension correction.}
#'   }
#'
#' @export
caps_forecast <- function(fit,
                          data,
                          horizons = c(7L, 14L, 21L, 28L),
                          alpha    = 0.05,
                          gamma_0  = 0.05,
                          method   = c("caps", "caps_static", "caps_marginal"),
                          n_sim    = 1000L,
                          ...) {

  method <- match.arg(method)

  if (!inherits(fit, "lfq_fit")) {
    cli::cli_abort("{.arg fit} must be an {.cls lfq_fit} object.")
  }
  if (!inherits(data, "lfq_data")) {
    cli::cli_abort("{.arg data} must be an {.cls lfq_data} object.")
  }

  lineages <- fit$lineages
  K        <- length(lineages)

  # ─── Stage 0: Generate parametric forecasts ──────────────────────────
  parametric_forecasts <- list()
  for (h in horizons) {
    fc <- tryCatch(
      forecast(fit, horizon = h, n_sim = n_sim),
      error = function(e) NULL
    )
    if (!is.null(fc)) {
      parametric_forecasts[[as.character(h)]] <- fc
    }
  }

  if (length(parametric_forecasts) == 0) {
    cli::cli_abort("Could not generate parametric forecasts for any horizon.")
  }

  # ─── Stage 1: Internal backtest for calibration ──────────────────────
  cli::cli_text("[CAPS] Running internal backtest for calibration...")

  bt <- tryCatch(
    backtest(data, engines = fit$engine, horizons = horizons,
             min_train = 42L, ...),
    error = function(e) {
      cli::cli_abort("Backtest failed: {e$message}")
    }
  )

  # bt is an lfq_backtest tibble directly
  # Columns: origin_date, target_date, horizon, engine, lineage,
  #          predicted, lower, upper, observed

  # Variance ratios R_hat(h)
  R_hat <- .caps_variance_ratios(bt, horizons)
  cli::cli_text("[CAPS] Variance ratios: {paste(sprintf('h=%dd R=%.2f', horizons, R_hat), collapse = ', ')}")

  # Residuals by horizon
  residuals_by_h <- .caps_residuals_by_horizon(bt, horizons)

  # ─── Stage 2: Horizon-specific CAPS parameters ──────────────────────
  caps_params <- list()
  for (i in seq_along(horizons)) {
    h     <- horizons[i]
    h_key <- as.character(h)
    resids <- residuals_by_h[[h_key]]

    if (is.null(resids) || length(resids) < 3L) {
      caps_params[[h_key]] <- list(
        horizon = h, R_hat = NA_real_, r_base = NA_real_,
        gamma_h = gamma_0, alpha_t = alpha, n_cal = 0L,
        q_base = NA_real_
      )
      next
    }

    Rh <- R_hat[i]
    if (is.na(Rh) || Rh <= 0) Rh <- 0.5

    # Conformal quantile
    n_cal   <- length(resids)
    q_level <- min(ceiling((1 - alpha) * (n_cal + 1)) / n_cal, 1)
    q_base  <- as.numeric(stats::quantile(abs(resids), probs = q_level,
                                           na.rm = TRUE))

    # Variance-ratio correction: widen when R < 1
    r_base_h <- q_base * sqrt(max(1, 1 / Rh))

    # Horizon-specific learning rate
    gamma_h <- if (method == "caps_static") 0 else gamma_0 * max(0, 1 - Rh)

    caps_params[[h_key]] <- list(
      horizon = h, R_hat = Rh, r_base = r_base_h,
      gamma_h = gamma_h, alpha_t = alpha,
      n_cal = n_cal, q_base = q_base
    )
  }

  # ─── Stage 3: Compositional dimension correction ────────────────────
  if (K > 2L && method == "caps") {
    psi_K <- (gamma((K - 1) / 2 + 1) / (pi^((K - 1) / 2)))^(1 / (K - 1))
  } else {
    psi_K <- 1
  }
  cli::cli_text("[CAPS] Dimension correction psi({K}) = {round(psi_K, 4)}")

  # ─── Stage 4: Apply CAPS intervals to parametric forecasts ──────────
  result_rows <- list()

  for (h_key in names(parametric_forecasts)) {
    h      <- as.integer(h_key)
    fc     <- parametric_forecasts[[h_key]]
    params <- caps_params[[h_key]]

    if (is.null(params) || is.na(params$r_base)) next

    # Extract forecast rows only (not fitted values)
    fc_future <- fc[fc$.type == "forecast", ]
    if (nrow(fc_future) == 0) next

    # CAPS radius
    r_caps <- params$r_base * psi_K

    # Build result tibble matching backtest column conventions
    for (ri in seq_len(nrow(fc_future))) {
      row <- fc_future[ri, ]
      result_rows[[length(result_rows) + 1L]] <- tibble::tibble(
        horizon       = h,
        lineage       = row$.lineage,
        target_date   = row$.date,
        predicted     = row$.median,
        lower_param   = row$.lower,
        upper_param   = row$.upper,
        caps_lower    = max(0, row$.median - r_caps),
        caps_upper    = min(1, row$.median + r_caps),
        caps_radius   = r_caps,
        caps_method   = method
      )
    }
  }

  forecasts_df <- if (length(result_rows) > 0) {
    dplyr::bind_rows(result_rows)
  } else {
    tibble::tibble()
  }

  # ─── Return ─────────────────────────────────────────────────────────
  structure(
    list(
      forecasts = forecasts_df,
      params    = caps_params,
      R_hat     = R_hat,
      psi_K     = psi_K,
      method    = method,
      alpha     = alpha,
      gamma_0   = gamma_0,
      horizons  = horizons,
      K         = K,
      lineages  = lineages
    ),
    class = "caps_forecast"
  )
}


#' Evaluate CAPS coverage on backtest data
#'
#' @param caps A \code{caps_forecast} object.
#' @param bt An \code{lfq_backtest} object to evaluate against.
#'
#' @return Data frame with coverage, width, and Winkler score by horizon,
#'   comparing CAPS, split conformal, and parametric intervals.
#' @export
evaluate_caps <- function(caps, bt) {

  if (!inherits(caps, "caps_forecast")) {
    cli::cli_abort("{.arg caps} must be a {.cls caps_forecast} object.")
  }
  if (!inherits(bt, "lfq_backtest")) {
    cli::cli_abort("{.arg bt} must be an {.cls lfq_backtest} object.")
  }

  # Split backtest into calibration (60%) and test (40%)
  origins     <- sort(unique(bt$origin_date))
  n_origins   <- length(origins)
  n_cal       <- floor(0.6 * n_origins)
  cal_origins <- origins[seq_len(n_cal)]
  test_origins <- origins[(n_cal + 1L):n_origins]

  cal_data  <- bt[bt$origin_date %in% cal_origins, ]
  test_data <- bt[bt$origin_date %in% test_origins, ]

  winkler_fn <- function(lo, hi, obs, a = caps$alpha) {
    w   <- hi - lo
    p_l <- ifelse(obs < lo, (2 / a) * (lo - obs), 0)
    p_u <- ifelse(obs > hi, (2 / a) * (obs - hi), 0)
    mean(w + p_l + p_u, na.rm = TRUE)
  }

  results <- list()

  for (h in caps$horizons) {
    h_key  <- as.character(h)
    test_h <- test_data[test_data$horizon == h, ]
    cal_h  <- cal_data[cal_data$horizon == h, ]

    if (nrow(test_h) < 3L || nrow(cal_h) < 3L) next

    # 1. Parametric
    para_cov   <- mean(test_h$observed >= test_h$lower &
                        test_h$observed <= test_h$upper, na.rm = TRUE)
    para_width <- mean(test_h$upper - test_h$lower, na.rm = TRUE)
    para_wink  <- winkler_fn(test_h$lower, test_h$upper, test_h$observed)

    # 2. Standard split conformal
    cal_resids <- abs(cal_h$predicted - cal_h$observed)
    nc         <- length(cal_resids)
    q_conf     <- as.numeric(stats::quantile(
      cal_resids, probs = min(1, ceiling((1 - caps$alpha) * (nc + 1)) / nc)
    ))
    conf_lo    <- pmax(0, test_h$predicted - q_conf)
    conf_hi    <- pmin(1, test_h$predicted + q_conf)
    conf_cov   <- mean(test_h$observed >= conf_lo &
                        test_h$observed <= conf_hi, na.rm = TRUE)
    conf_width <- mean(conf_hi - conf_lo, na.rm = TRUE)
    conf_wink  <- winkler_fn(conf_lo, conf_hi, test_h$observed)

    # 3. CAPS
    params <- caps$params[[h_key]]
    if (is.null(params) || is.na(params$r_base)) next
    r_caps  <- params$r_base * caps$psi_K
    caps_lo <- pmax(0, test_h$predicted - r_caps)
    caps_hi <- pmin(1, test_h$predicted + r_caps)
    caps_cov   <- mean(test_h$observed >= caps_lo &
                        test_h$observed <= caps_hi, na.rm = TRUE)
    caps_width <- mean(caps_hi - caps_lo, na.rm = TRUE)
    caps_wink  <- winkler_fn(caps_lo, caps_hi, test_h$observed)

    results[[h_key]] <- tibble::tibble(
      horizon       = h,
      R_hat         = params$R_hat,
      para_cov      = para_cov,
      para_width    = para_width,
      para_winkler  = para_wink,
      conf_cov      = conf_cov,
      conf_width    = conf_width,
      conf_winkler  = conf_wink,
      caps_cov      = caps_cov,
      caps_width    = caps_width,
      caps_winkler  = caps_wink,
      n_test        = nrow(test_h)
    )
  }

  dplyr::bind_rows(results)
}


#' @export
print.caps_forecast <- function(x, ...) {
  cli::cli_h3("CAPS forecast")
  cli::cli_text("Method: {x$method} | Alpha: {x$alpha} | K: {x$K} lineages")
  cli::cli_text("Dimension correction psi({x$K}) = {round(x$psi_K, 4)}")
  cli::cli_text("")
  for (h_key in names(x$params)) {
    p <- x$params[[h_key]]
    if (!is.na(p$R_hat)) {
      cli::cli_text("  h={h_key}d: R_hat={round(p$R_hat, 3)}, r_base={round(p$r_base, 4)}, gamma={round(p$gamma_h, 4)}, n_cal={p$n_cal}")
    }
  }
  if (nrow(x$forecasts) > 0) {
    cli::cli_text("")
    cli::cli_text("{nrow(x$forecasts)} forecast rows across {length(unique(x$forecasts$horizon))} horizons")
  }
  invisible(x)
}


# ── Internal helpers ───────────────────────────────────────────────────

#' Compute variance ratios by horizon from backtest tibble
#' @noRd
.caps_variance_ratios <- function(bt, horizons) {
  R_hat <- numeric(length(horizons))
  names(R_hat) <- as.character(horizons)

  for (i in seq_along(horizons)) {
    h    <- horizons[i]
    bt_h <- bt[bt$horizon == h &
                 !is.na(bt$predicted) & !is.na(bt$observed) &
                 !is.na(bt$lower) & !is.na(bt$upper), ]

    if (nrow(bt_h) < 5L) {
      R_hat[i] <- NA_real_
      next
    }

    pred_var <- mean(((bt_h$upper - bt_h$lower) / (2 * stats::qnorm(0.975)))^2,
                     na.rm = TRUE)
    obs_mse  <- mean((bt_h$predicted - bt_h$observed)^2, na.rm = TRUE)

    R_hat[i] <- if (obs_mse > 1e-15) pred_var / obs_mse else NA_real_
  }

  R_hat
}

#' Compute residuals by horizon from backtest tibble
#' @noRd
.caps_residuals_by_horizon <- function(bt, horizons) {
  res <- list()
  for (h in horizons) {
    bt_h <- bt[bt$horizon == h &
                 !is.na(bt$predicted) & !is.na(bt$observed), ]
    if (nrow(bt_h) > 0) {
      res[[as.character(h)]] <- bt_h$predicted - bt_h$observed
    }
  }
  res
}
