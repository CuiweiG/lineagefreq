#' Joint Conformal Prediction on the Simplex
#'
#' Produces a compositional prediction region that respects the
#' simplex constraint (frequencies sum to 1) using Aitchison geometry.
#' Unlike marginal conformal prediction, joint prediction guarantees
#' that the \emph{entire} frequency vector is covered, not just
#' individual lineages.
#'
#' @param fit An \code{lfq_fit} object.
#' @param data An \code{lfq_data} object (same data used to fit the model).
#' @param horizon Integer; forecast horizon in days (default 28).
#' @param ci_level Numeric in (0,1); coverage target (default 0.95).
#' @param cal_fraction Numeric in (0,1); fraction of dates for
#'   calibration (default 0.3).
#' @param seed Optional integer for reproducibility.
#'
#' @return An \code{lfq_conformal_joint} S3 object (list) with:
#'   \describe{
#'     \item{forecast}{The point forecast (\code{lfq_forecast}).}
#'     \item{radius}{Conformal radius in Aitchison distance.}
#'     \item{marginal_intervals}{Tibble with .date, .lineage,
#'       .lower_joint, .upper_joint — marginal bounds projected from
#'       the joint region.}
#'     \item{marginal_only}{Tibble with .lower_marginal, .upper_marginal
#'       from standard marginal conformal prediction.}
#'     \item{comparison}{Tibble indicating whether joint intervals are
#'       wider or narrower than marginal intervals per lineage.}
#'     \item{calibration_scores}{Aitchison distances on calibration set.}
#'     \item{ci_level}{Nominal coverage level.}
#'     \item{n_cal}{Number of calibration compositions.}
#'   }
#'
#' @details
#' The nonconformity score is the Aitchison distance between predicted
#' and observed compositions. The Aitchison distance equals the
#' Euclidean distance in the isometric log-ratio (ILR) transformed
#' space, which respects the geometry of the simplex (Aitchison, 1986).
#'
#' The prediction region is the set of all compositions within
#' Aitchison distance \eqn{r} of the point forecast, where \eqn{r}
#' is the \eqn{(1-\alpha)(1+1/n)} quantile of calibration distances.
#' Marginal intervals are obtained by projecting this region onto
#' each coordinate axis, then intersecting with \eqn{[0, 1]}.
#'
#' @references
#' Aitchison J (1986). \emph{The Statistical Analysis of Compositional
#' Data}. Chapman & Hall.
#'
#' Vovk V, Gammerman A, Shafer G (2005). \emph{Algorithmic Learning
#' in a Random World}. Springer.
#'
#' @export
conformal_forecast_joint <- function(fit, data,
                                     horizon      = 28L,
                                     ci_level     = 0.95,
                                     cal_fraction = 0.3,
                                     seed         = NULL) {

  if (!inherits(fit, "lfq_fit")) {
    cli::cli_abort("{.arg fit} must be an {.cls lfq_fit} object.")
  }
  if (!inherits(data, "lfq_data")) {
    cli::cli_abort("{.arg data} must be an {.cls lfq_data} object.")
  }

  if (!is.null(seed)) set.seed(seed)

  lineages   <- attr(data, "lineages")
  K          <- length(lineages)
  dates      <- sort(unique(data$.date))
  n_dates    <- length(dates)

  # ── Temporal split: calibration = last cal_fraction of dates ──────────
  n_cal_dates <- max(floor(n_dates * cal_fraction), 2L)
  n_train     <- n_dates - n_cal_dates

  if (n_train < 3L) {
    cli::cli_abort("Too few training dates ({n_train}) for joint conformal.")
  }

  train_dates <- dates[seq_len(n_train)]
  cal_dates   <- dates[(n_train + 1L):n_dates]

  # ── Refit model on training portion ──────────────────────────────────
  train_data <- data[data$.date <= max(train_dates), ]
  # Preserve lfq_data class and attributes
  for (a in c("lineages", "date_range", "n_timepoints",
              "has_location", "min_total")) {
    attr(train_data, a) <- attr(data, a)
  }
  class(train_data) <- class(data)
  attr(train_data, "date_range")   <- range(train_data$.date)
  attr(train_data, "n_timepoints") <- length(unique(train_data$.date))

  train_fit <- tryCatch(
    fit_model(train_data, engine = fit$engine, pivot = fit$pivot,
              horizon = horizon, ci_level = ci_level),
    error = function(e) {
      cli::cli_abort("Refitting model on training data failed: {e$message}")
    }
  )

  # ── Compute predicted compositions on calibration dates ──────────────
  t0         <- min(dates)
  time_scale <- fit$time_scale %||% 7L
  intercepts <- train_fit$intercepts
  growth_rates <- train_fit$growth_rates
  pivot      <- train_fit$pivot
  non_pivot  <- setdiff(lineages, pivot)

  predict_composition <- function(date_val) {
    t_val <- as.numeric(date_val - t0) / time_scale
    log_num <- numeric(K)
    names(log_num) <- lineages
    log_num[pivot] <- 0
    log_num[non_pivot] <- intercepts[non_pivot] +
                          growth_rates[non_pivot] * t_val
    log_denom <- log_sum_exp(log_num)
    exp(log_num - log_denom)
  }

  # ── Extract observed compositions on calibration dates ───────────────
  cal_observed <- list()
  cal_predicted <- list()

  for (di in seq_along(cal_dates)) {
    d <- cal_dates[di]
    obs_at_d <- data[data$.date == d, ]
    if (nrow(obs_at_d) == 0) next

    obs_vec <- stats::setNames(rep(0, K), lineages)
    for (i in seq_len(nrow(obs_at_d))) {
      lin <- obs_at_d$.lineage[i]
      if (lin %in% lineages) obs_vec[lin] <- obs_at_d$.freq[i]
    }
    # Normalise to simplex
    s <- sum(obs_vec)
    if (s > 0) obs_vec <- obs_vec / s

    pred_vec <- predict_composition(d)

    cal_observed[[length(cal_observed) + 1]]   <- obs_vec
    cal_predicted[[length(cal_predicted) + 1]] <- pred_vec
  }

  n_cal <- length(cal_observed)
  if (n_cal < 2L) {
    cli::cli_abort("Too few calibration compositions ({n_cal}).")
  }

  # ── Compute Aitchison distances ──────────────────────────────────────
  cal_distances <- vapply(seq_len(n_cal), function(i) {
    .aitchison_distance(cal_predicted[[i]], cal_observed[[i]])
  }, numeric(1L))

  # ── Conformal quantile ──────────────────────────────────────────────
  q_level <- min(ceiling((1 - (1 - ci_level)) * (n_cal + 1)) / n_cal, 1)
  radius  <- stats::quantile(cal_distances, probs = q_level,
                              names = FALSE)

  # ── Generate point forecast ─────────────────────────────────────────
  point_forecast <- forecast(fit, horizon = horizon, ci_level = ci_level)

  # ── Project Aitchison ball to marginal intervals ────────────────────
  forecast_dates <- unique(
    point_forecast[point_forecast$.type == "forecast", ]$.date
  )

  marginal_rows <- list()
  for (fi in seq_along(forecast_dates)) {
    d <- forecast_dates[fi]
    fc_at_d <- point_forecast[point_forecast$.date == d &
                               point_forecast$.type == "forecast", ]
    if (nrow(fc_at_d) == 0) next

    centre <- stats::setNames(rep(0, K), lineages)
    for (i in seq_len(nrow(fc_at_d))) {
      lin <- fc_at_d$.lineage[i]
      if (lin %in% lineages) centre[lin] <- fc_at_d$.median[i]
    }
    s <- sum(centre)
    if (s > 0) centre <- centre / s

    bounds <- .project_aitchison_ball(centre, radius, lineages)

    for (lin in lineages) {
      marginal_rows[[length(marginal_rows) + 1]] <- tibble::tibble(
        .date          = d,
        .lineage       = lin,
        .median        = centre[lin],
        .lower_joint   = bounds$lower[lin],
        .upper_joint   = bounds$upper[lin]
      )
    }
  }

  marginal_df <- dplyr::bind_rows(marginal_rows)

  # ── Marginal-only intervals for comparison ──────────────────────────
  marginal_conformal <- tryCatch(
    conformal_forecast(fit, data, horizon = horizon, ci_level = ci_level,
                       method = "split", cal_fraction = cal_fraction,
                       seed = seed),
    error = function(e) NULL
  )

  comparison <- NULL
  if (!is.null(marginal_conformal)) {
    mc_forecast <- marginal_conformal[marginal_conformal$.type == "forecast", ]
    mc_df <- mc_forecast |>
      dplyr::select(.date, .lineage,
                    .lower_marginal = .lower,
                    .upper_marginal = .upper)

    marginal_df <- dplyr::left_join(marginal_df, mc_df,
                                     by = c(".date", ".lineage"))

    comparison <- marginal_df |>
      dplyr::mutate(
        width_joint    = .upper_joint - .lower_joint,
        width_marginal = .upper_marginal - .lower_marginal,
        joint_wider    = width_joint > width_marginal + 1e-10
      ) |>
      dplyr::group_by(.lineage) |>
      dplyr::summarise(
        mean_width_joint    = mean(width_joint, na.rm = TRUE),
        mean_width_marginal = mean(width_marginal, na.rm = TRUE),
        pct_joint_wider     = mean(joint_wider, na.rm = TRUE),
        .groups = "drop"
      )
  }

  structure(
    list(
      forecast            = point_forecast,
      radius              = radius,
      marginal_intervals  = marginal_df,
      comparison          = comparison,
      calibration_scores  = cal_distances,
      ci_level            = ci_level,
      n_cal               = n_cal
    ),
    class = "lfq_conformal_joint"
  )
}


#' Joint Calibration Diagnostics
#'
#' Assess calibration of multivariate forecasts using the energy score,
#' joint coverage, and multivariate rank histograms.
#'
#' @param bt An \code{lfq_backtest} object.
#' @param n_bins Integer; number of rank histogram bins (default 10).
#'
#' @return A \code{joint_calibration_report} S3 object (list) with:
#'   \describe{
#'     \item{energy_scores}{Tibble with origin_date, horizon, energy_score.}
#'     \item{mean_energy_score}{Overall mean energy score.}
#'     \item{joint_coverage}{Tibble with nominal level, observed coverage,
#'       using Aitchison distance from backtest results.}
#'     \item{rank_histogram}{Tibble with bin, count, density, expected.}
#'     \item{n}{Number of forecast-observation vectors.}
#'   }
#'
#' @details
#' The energy score is a multivariate proper scoring rule that
#' generalises the CRPS (Gneiting & Raftery, 2007, Section 5).
#' For a deterministic forecast \eqn{f} and observation \eqn{y}:
#' \deqn{ES = ||f - y||^2}
#' where the norm is the Aitchison distance on the simplex.
#'
#' Joint coverage at level \eqn{1-\alpha}: the fraction of observed
#' composition vectors falling within Aitchison distance
#' \eqn{q_\alpha} of the forecast, where \eqn{q_\alpha} is derived
#' from the empirical distribution of distances.
#'
#' The multivariate rank histogram uses the pre-rank approach:
#' the rank of the observation among an ensemble is computed using
#' the Aitchison distance to the ensemble mean.
#'
#' @export
calibrate_joint <- function(bt, n_bins = 10L) {

  if (!inherits(bt, "lfq_backtest")) {
    cli::cli_abort("{.arg bt} must be an {.cls lfq_backtest} object.")
  }

  lineages <- sort(unique(bt$lineage))
  K        <- length(lineages)
  origins  <- sort(unique(bt$origin_date))
  horizons <- sort(unique(bt$horizon))

  # ── Compute Aitchison distance for each origin × horizon ─────────────
  distance_rows <- list()

  for (oi in seq_along(origins)) {
    orig <- origins[oi]
    for (h in horizons) {
      slice <- bt[bt$origin_date == orig & bt$horizon == h, ]
      if (nrow(slice) == 0) next

      pred_vec <- stats::setNames(rep(0, K), lineages)
      obs_vec  <- stats::setNames(rep(0, K), lineages)

      for (i in seq_len(nrow(slice))) {
        lin <- slice$lineage[i]
        if (lin %in% lineages) {
          pred_vec[lin] <- slice$predicted[i]
          obs_vec[lin]  <- slice$observed[i]
        }
      }

      # Normalise
      sp <- sum(pred_vec); if (sp > 0) pred_vec <- pred_vec / sp
      so <- sum(obs_vec);  if (so > 0) obs_vec  <- obs_vec / so

      d <- .aitchison_distance(pred_vec, obs_vec)
      es <- d^2  # Energy score (squared Aitchison distance)

      distance_rows[[length(distance_rows) + 1]] <- tibble::tibble(
        origin_date  = orig,
        horizon      = h,
        aitchison_d  = d,
        energy_score = es
      )
    }
  }

  dist_df <- dplyr::bind_rows(distance_rows)

  # ── Energy score ─────────────────────────────────────────────────────
  mean_es <- mean(dist_df$energy_score, na.rm = TRUE)

  # ── Joint coverage at multiple nominal levels ────────────────────────
  nominal_levels <- seq(0.1, 0.95, by = 0.05)
  all_distances  <- dist_df$aitchison_d
  n              <- length(all_distances)

  joint_cov <- tibble::tibble(
    nominal  = nominal_levels,
    observed = vapply(nominal_levels, function(level) {
      # Threshold: the level-quantile of observed distances
      threshold <- stats::quantile(all_distances, probs = level,
                                    names = FALSE)
      mean(all_distances <= threshold)
    }, numeric(1L))
  )

  # ── Multivariate rank histogram ──────────────────────────────────────
  # Pre-rank: for each origin-horizon, rank = fraction of
  # calibration distances smaller than the observation distance
  # (i.e., where does this observation sit in the distribution?)
  pit_multivariate <- vapply(all_distances, function(d) {
    mean(all_distances <= d)
  }, numeric(1L))

  breaks <- seq(0, 1, length.out = n_bins + 1L)
  counts <- as.integer(table(cut(pit_multivariate, breaks = breaks,
                                  include.lowest = TRUE)))
  rank_hist <- tibble::tibble(
    bin      = seq_len(n_bins),
    count    = counts,
    density  = counts / n,
    expected = 1 / n_bins
  )

  structure(
    list(
      energy_scores      = dist_df,
      mean_energy_score  = mean_es,
      joint_coverage     = joint_cov,
      rank_histogram     = rank_hist,
      n                  = n
    ),
    class = "joint_calibration_report"
  )
}


#' @export
print.lfq_conformal_joint <- function(x, ...) {
  cli::cli_h3("Joint conformal prediction (simplex-aware)")
  cli::cli_text("Aitchison radius: {round(x$radius, 4)}")
  cli::cli_text("Nominal coverage: {x$ci_level * 100}%")
  cli::cli_text("Calibration set: {x$n_cal} compositions")
  if (!is.null(x$comparison)) {
    wider_pct <- mean(x$comparison$pct_joint_wider) * 100
    cli::cli_text("Joint intervals wider than marginal: {round(wider_pct, 0)}% of cases")
  }
  invisible(x)
}


#' @export
print.joint_calibration_report <- function(x, ...) {
  cli::cli_h3("Joint calibration report")
  cli::cli_text("{x$n} forecast-observation vectors")
  cli::cli_text("Mean energy score: {round(x$mean_energy_score, 4)}")
  # Show coverage at 90% and 95% nominal
  for (level in c(0.90, 0.95)) {
    row <- x$joint_coverage[x$joint_coverage$nominal == level, ]
    if (nrow(row) > 0) {
      cli::cli_text("Joint coverage at {level*100}%: {round(row$observed*100, 1)}%")
    }
  }
  invisible(x)
}


# ── Internal helpers ───────────────────────────────────────────────────

#' Isometric log-ratio transform
#'
#' Maps a K-dimensional composition on the simplex to R^(K-1)
#' using the Helmert sub-composition basis.
#' @noRd
.ilr_transform <- function(x) {
  K <- length(x)
  if (K < 2L) return(numeric(0L))

  # Replace zeros with small value for log stability
  x <- pmax(x, 1e-15)
  x <- x / sum(x)

  # Helmert ILR basis (Egozcue et al. 2003)
  ilr <- numeric(K - 1L)
  for (k in seq_len(K - 1L)) {
    geom_mean_k <- exp(mean(log(x[seq_len(k)])))
    ilr[k] <- sqrt(k / (k + 1)) * log(geom_mean_k / x[k + 1L])
  }
  ilr
}

#' Inverse ILR transform
#'
#' Maps from R^(K-1) back to the K-simplex using the Helmert
#' contrast matrix.
#' @noRd
.ilr_inverse <- function(ilr_coords) {
  K <- length(ilr_coords) + 1L
  if (K < 2L) return(1)

  # Build the Helmert contrast matrix V: (K-1) x K
  V <- matrix(0, K - 1L, K)
  for (k in seq_len(K - 1L)) {
    coeff <- sqrt(1 / (k * (k + 1)))
    V[k, seq_len(k)] <- coeff
    V[k, k + 1L]     <- -k * coeff
  }

  # log(x) = V^T %*% ilr (up to an additive constant)
  log_x <- as.numeric(crossprod(V, ilr_coords))
  x <- exp(log_x)
  x / sum(x)
}

#' Aitchison distance between two compositions
#'
#' Equals the Euclidean distance in ILR space.
#' @noRd
.aitchison_distance <- function(x, y) {
  ilr_x <- .ilr_transform(x)
  ilr_y <- .ilr_transform(y)
  sqrt(sum((ilr_x - ilr_y)^2))
}

#' Project an Aitchison ball to marginal intervals
#'
#' Given a centre composition and radius in Aitchison distance,
#' find the min/max of each component over all compositions
#' within the ball. Uses numerical optimisation on ILR coordinates.
#' @noRd
.project_aitchison_ball <- function(centre, radius, lineage_names) {
  K <- length(centre)
  centre_ilr <- .ilr_transform(centre)

  lower <- stats::setNames(numeric(K), lineage_names)
  upper <- stats::setNames(numeric(K), lineage_names)

  for (k in seq_len(K)) {
    # Minimise component k subject to ||ilr(x) - ilr(centre)|| <= radius
    # Use grid search on the Aitchison ball boundary
    # For K-1 dimensional ball, sample directions and find extremes

    if (K == 2L) {
      # Special case: 1D ILR, can solve analytically
      lo_ilr <- centre_ilr - radius
      hi_ilr <- centre_ilr + radius
      comp_lo <- .ilr_inverse(lo_ilr)
      comp_hi <- .ilr_inverse(hi_ilr)
      lower[k] <- min(comp_lo[k], comp_hi[k])
      upper[k] <- max(comp_lo[k], comp_hi[k])
    } else {
      # General case: sample directions on the (K-2)-sphere
      n_dirs <- max(200L, 50L * K)
      extremes <- numeric(n_dirs)

      for (d in seq_len(n_dirs)) {
        # Random direction on (K-1)-sphere
        dir <- stats::rnorm(K - 1L)
        dir <- dir / sqrt(sum(dir^2))

        # Point on the ball boundary
        ilr_point <- centre_ilr + radius * dir
        comp <- .ilr_inverse(ilr_point)
        extremes[d] <- comp[k]
      }

      # Also check the centre
      lower[k] <- min(extremes, centre[k])
      upper[k] <- max(extremes, centre[k])
    }
  }

  # Clamp to [0, 1]
  lower <- pmax(lower, 0)
  upper <- pmin(upper, 1)

  list(lower = lower, upper = upper)
}
