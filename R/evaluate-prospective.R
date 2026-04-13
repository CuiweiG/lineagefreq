#' Pseudo-Prospective Evaluation of Conformal Prediction
#'
#' Simulates real-time deployment: at each origin, all preceding
#' origins form the conformal calibration set, and the forecast at
#' the next origin is evaluated. No future data is ever used for
#' calibration — this is a genuine expanding-window evaluation.
#'
#' @param data An \code{lfq_data} object.
#' @param engine Character; estimation engine (default \code{"mlr"}).
#' @param horizons Integer vector; forecast horizons in days
#'   (default \code{c(7L, 14L)}).
#' @param min_train Integer; minimum training window in days
#'   (default 42).
#' @param min_cal Integer; minimum calibration origins before
#'   conformal intervals are computed (default 5).
#' @param ci_level Numeric in (0,1); nominal coverage (default 0.95).
#' @param gamma Numeric; ACI learning rate (default 0.05).
#' @param ... Passed to \code{fit_model()}.
#'
#' @return An \code{lfq_prospective} S3 object (list) with:
#'   \describe{
#'     \item{results}{Tibble with origin_date, target_date, horizon,
#'       lineage, predicted, observed, lower_param, upper_param,
#'       lower_static, upper_static, lower_aci, upper_aci,
#'       radius_static, radius_aci, alpha_aci.}
#'     \item{coverage}{Tibble summarising cumulative coverage by
#'       method (parametric, static conformal, ACI) at each step.}
#'     \item{summary}{Tibble with method, overall coverage, mean
#'       width, Winkler score.}
#'     \item{n_origins}{Total number of evaluation origins.}
#'   }
#'
#' @details
#' At each origin \eqn{t_i} with \eqn{i \geq} \code{min_cal}:
#' \enumerate{
#'   \item Fit the model on data up to \eqn{t_i}.
#'   \item Compute residuals at all previous origins
#'         \eqn{t_1, \ldots, t_{i-1}} (the calibration set).
#'   \item Static conformal: radius = \eqn{(1-\alpha)(1+1/n)}
#'         quantile of absolute calibration residuals.
#'   \item ACI: radius adjusted by online update of \eqn{\alpha_t}.
#'   \item Evaluate coverage at \eqn{t_{i+1}}.
#' }
#'
#' This differs from \code{conformal_forecast()} which uses a
#' fixed temporal split. Here, the calibration set expands over time,
#' mimicking a surveillance agency accumulating validation data.
#'
#' @export
evaluate_prospective <- function(data,
                                  engine    = "mlr",
                                  horizons  = c(7L, 14L),
                                  min_train = 42L,
                                  min_cal   = 5L,
                                  ci_level  = 0.95,
                                  gamma     = 0.05,
                                  ...) {

  if (!inherits(data, "lfq_data")) {
    cli::cli_abort("{.arg data} must be an {.cls lfq_data} object.")
  }

  alpha     <- 1 - ci_level
  lineages  <- attr(data, "lineages")
  dates     <- sort(unique(data$.date))
  n_dates   <- length(dates)
  min_date  <- min(dates)
  max_date  <- max(dates)

  # ── Determine all valid origin dates ────────────────────────────────
  eligible <- dates[as.numeric(dates - min_date) >= min_train]
  # Exclude last date (no observation to evaluate against)
  eligible <- eligible[eligible < max_date]

  if (length(eligible) < min_cal + 1L) {
    cli::cli_abort("Too few eligible origins ({length(eligible)}). Need at least {min_cal + 1}.")
  }

  cli::cli_text("Prospective evaluation: {length(eligible)} origins, {length(horizons)} horizons")

  # ── Storage ─────────────────────────────────────────────────────────
  all_results   <- list()
  cal_residuals <- list()  # list of lists, one per horizon
  for (h in horizons) cal_residuals[[as.character(h)]] <- numeric(0)

  alpha_aci <- alpha  # ACI running miscoverage rate

  # ── Main loop: expanding window ─────────────────────────────────────
  for (idx in seq_along(eligible)) {
    origin <- eligible[idx]

    # Fit model on data up to origin
    train_data <- data[data$.date <= origin, ]
    for (a in c("lineages", "date_range", "n_timepoints",
                "has_location", "min_total")) {
      attr(train_data, a) <- attr(data, a)
    }
    class(train_data) <- class(data)
    attr(train_data, "date_range")   <- range(train_data$.date)
    attr(train_data, "n_timepoints") <- length(unique(train_data$.date))

    fit <- tryCatch(
      fit_model(train_data, engine = engine, ...),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    # Generate parametric forecast
    fc <- tryCatch(
      forecast(fit, horizon = max(horizons), ci_level = ci_level),
      error = function(e) NULL
    )
    if (is.null(fc)) next

    fc_future <- fc[fc$.type == "forecast", ]

    for (h in horizons) {
      h_key      <- as.character(h)
      target_date <- origin + h

      # Find observed frequencies at target
      obs_at_target <- data[data$.date == target_date, ]
      if (nrow(obs_at_target) == 0) {
        # Try closest date within 3 days
        close_dates <- dates[abs(as.numeric(dates - target_date)) <= 3]
        if (length(close_dates) > 0) {
          obs_at_target <- data[data$.date == close_dates[1], ]
          target_date <- close_dates[1]
        } else {
          next
        }
      }

      # Find forecast closest to target_date
      fc_dates <- unique(fc_future$.date)
      if (length(fc_dates) == 0) next
      closest_fc_date <- fc_dates[which.min(abs(as.numeric(fc_dates - target_date)))]
      fc_at_h <- fc_future[fc_future$.date == closest_fc_date, ]

      for (lin in lineages) {
        fc_row  <- fc_at_h[fc_at_h$.lineage == lin, ]
        obs_row <- obs_at_target[obs_at_target$.lineage == lin, ]

        if (nrow(fc_row) == 0 || nrow(obs_row) == 0) next

        pred <- fc_row$.median[1]
        obs  <- obs_row$.freq[1]
        lo_p <- fc_row$.lower[1]
        hi_p <- fc_row$.upper[1]

        # Compute residual for this origin
        resid <- abs(pred - obs)

        # ── Static conformal ──────────────────────────────────────
        cal_res <- cal_residuals[[h_key]]
        n_cal   <- length(cal_res)

        if (n_cal >= min_cal) {
          q_level <- min(ceiling((1 - alpha) * (n_cal + 1)) / n_cal, 1)
          radius_static <- stats::quantile(cal_res, probs = q_level,
                                            names = FALSE)
          lo_s <- pmax(pred - radius_static, 0)
          hi_s <- pmin(pred + radius_static, 1)
        } else {
          radius_static <- NA_real_
          lo_s <- NA_real_
          hi_s <- NA_real_
        }

        # ── ACI conformal ─────────────────────────────────────────
        if (n_cal >= min_cal) {
          q_aci <- min(ceiling((1 - alpha_aci) * (n_cal + 1)) / n_cal, 1)
          q_aci <- max(q_aci, 0.01)
          radius_aci <- stats::quantile(cal_res, probs = q_aci,
                                         names = FALSE)
          lo_a <- pmax(pred - radius_aci, 0)
          hi_a <- pmin(pred + radius_aci, 1)
        } else {
          radius_aci <- NA_real_
          lo_a <- NA_real_
          hi_a <- NA_real_
        }

        all_results[[length(all_results) + 1]] <- tibble::tibble(
          origin_date   = origin,
          target_date   = target_date,
          horizon       = h,
          lineage       = lin,
          predicted     = pred,
          observed      = obs,
          lower_param   = lo_p,
          upper_param   = hi_p,
          lower_static  = lo_s,
          upper_static  = hi_s,
          lower_aci     = lo_a,
          upper_aci     = hi_a,
          radius_static = radius_static,
          radius_aci    = radius_aci,
          alpha_aci     = alpha_aci,
          n_cal         = n_cal
        )

        # Update calibration residuals (add this origin's residual)
        cal_residuals[[h_key]] <- c(cal_res, resid)
      }

      # ── ACI online update ─────────────────────────────────────────
      # Update alpha based on whether the observation was covered by ACI
      if (n_cal >= min_cal) {
        covered_aci <- all(
          obs_at_target$.freq[obs_at_target$.lineage %in% lineages] >=
            (fc_at_h$.median - radius_aci)[fc_at_h$.lineage %in% lineages] &
          obs_at_target$.freq[obs_at_target$.lineage %in% lineages] <=
            (fc_at_h$.median + radius_aci)[fc_at_h$.lineage %in% lineages],
          na.rm = TRUE
        )
        # Gibbs & Candès update: alpha_t+1 = alpha_t + gamma*(alpha - err_t)
        err_t <- if (covered_aci) 0 else 1
        alpha_aci <- alpha_aci + gamma * (alpha - err_t)
        alpha_aci <- max(0.001, min(alpha_aci, 0.5))
      }
    }
  }

  results_df <- dplyr::bind_rows(all_results)

  if (nrow(results_df) == 0) {
    cli::cli_abort("No valid forecast-observation pairs produced.")
  }

  # ── Compute cumulative coverage ───────────────────────────────────
  results_df <- results_df |>
    dplyr::mutate(
      covered_param  = observed >= lower_param & observed <= upper_param,
      covered_static = !is.na(lower_static) &
                       observed >= lower_static & observed <= upper_static,
      covered_aci    = !is.na(lower_aci) &
                       observed >= lower_aci & observed <= upper_aci,
      width_param    = upper_param - lower_param,
      width_static   = upper_static - lower_static,
      width_aci      = upper_aci - lower_aci
    )

  # Cumulative coverage over origins (averaging across lineages per origin)
  origin_cov <- results_df |>
    dplyr::filter(!is.na(lower_static)) |>
    dplyr::group_by(origin_date, horizon) |>
    dplyr::summarise(
      cov_param  = mean(covered_param, na.rm = TRUE),
      cov_static = mean(covered_static, na.rm = TRUE),
      cov_aci    = mean(covered_aci, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(origin_date) |>
    dplyr::mutate(
      cum_param  = cumsum(cov_param) / seq_len(dplyr::n()),
      cum_static = cumsum(cov_static) / seq_len(dplyr::n()),
      cum_aci    = cumsum(cov_aci) / seq_len(dplyr::n()),
      step       = seq_len(dplyr::n())
    )

  # ── Summary statistics ───────────────────────────────────────────
  valid <- results_df |> dplyr::filter(!is.na(lower_static))
  winkler <- function(lo, hi, obs, alpha) {
    w <- hi - lo
    below <- pmax(0, lo - obs)
    above <- pmax(0, obs - hi)
    mean(w + (2 / alpha) * (below + above), na.rm = TRUE)
  }

  summary_df <- tibble::tibble(
    method = c("Parametric", "Static conformal", "ACI"),
    coverage = c(
      mean(valid$covered_param, na.rm = TRUE),
      mean(valid$covered_static, na.rm = TRUE),
      mean(valid$covered_aci, na.rm = TRUE)
    ),
    mean_width = c(
      mean(valid$width_param, na.rm = TRUE),
      mean(valid$width_static, na.rm = TRUE),
      mean(valid$width_aci, na.rm = TRUE)
    ),
    winkler_score = c(
      winkler(valid$lower_param, valid$upper_param, valid$observed, alpha),
      winkler(valid$lower_static, valid$upper_static, valid$observed, alpha),
      winkler(valid$lower_aci, valid$upper_aci, valid$observed, alpha)
    )
  )

  structure(
    list(
      results    = results_df,
      coverage   = origin_cov,
      summary    = summary_df,
      n_origins  = length(unique(results_df$origin_date)),
      ci_level   = ci_level,
      gamma      = gamma
    ),
    class = "lfq_prospective"
  )
}


#' @export
print.lfq_prospective <- function(x, ...) {
  cli::cli_h3("Prospective conformal evaluation")
  cli::cli_text("{x$n_origins} origins, {nrow(x$results)} forecast-obs pairs")
  cli::cli_text("Nominal coverage: {x$ci_level * 100}%")
  cli::cli_text("")
  for (i in seq_len(nrow(x$summary))) {
    s <- x$summary[i, ]
    cli::cli_text("{s$method}: coverage={round(s$coverage*100,1)}%, width={round(s$mean_width,4)}, Winkler={round(s$winkler_score,4)}")
  }
  invisible(x)
}


#' Plot Prospective Evaluation Results
#'
#' @param x An \code{lfq_prospective} object.
#' @param type Character; one of \code{"coverage"} (cumulative
#'   coverage over time), \code{"radius"} (conformal radius over
#'   time), or \code{"comparison"} (retrospective vs prospective).
#' @param ... Ignored.
#'
#' @return A ggplot object.
#' @export
plot.lfq_prospective <- function(x, type = c("coverage", "radius",
                                              "comparison"), ...) {
  type <- match.arg(type)

  if (type == "coverage") {
    cov_long <- x$coverage |>
      tidyr::pivot_longer(
        cols = c(cum_param, cum_static, cum_aci),
        names_to = "method", values_to = "coverage"
      ) |>
      dplyr::mutate(
        method = dplyr::case_when(
          method == "cum_param"  ~ "Parametric",
          method == "cum_static" ~ "Static conformal",
          method == "cum_aci"    ~ "ACI"
        )
      )

    ggplot2::ggplot(cov_long,
                    ggplot2::aes(x = step, y = coverage, colour = method)) +
      ggplot2::geom_line(linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = x$ci_level, linetype = "dashed",
                          colour = "grey50") +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::labs(x = "Evaluation step", y = "Cumulative coverage",
                    colour = NULL) +
      ggplot2::theme_classic()

  } else if (type == "radius") {
    radius_df <- x$results |>
      dplyr::filter(!is.na(radius_static)) |>
      dplyr::group_by(origin_date) |>
      dplyr::summarise(
        static = mean(radius_static, na.rm = TRUE),
        aci    = mean(radius_aci, na.rm = TRUE),
        .groups = "drop"
      ) |>
      tidyr::pivot_longer(cols = c(static, aci),
                          names_to = "method", values_to = "radius") |>
      dplyr::mutate(
        method = dplyr::case_when(
          method == "static" ~ "Static conformal",
          method == "aci"    ~ "ACI"
        )
      )

    ggplot2::ggplot(radius_df,
                    ggplot2::aes(x = origin_date, y = radius, colour = method)) +
      ggplot2::geom_line(linewidth = 0.6) +
      ggplot2::labs(x = "Origin date", y = "Conformal radius",
                    colour = NULL) +
      ggplot2::theme_classic()

  } else {
    # comparison: bar chart of coverage by method
    ggplot2::ggplot(x$summary,
                    ggplot2::aes(x = method, y = coverage,
                                 fill = method)) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::geom_hline(yintercept = x$ci_level, linetype = "dashed") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(),
                                  limits = c(0, 1)) +
      ggplot2::labs(x = NULL, y = "Coverage", fill = NULL) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none")
  }
}
