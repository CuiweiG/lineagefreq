#' Rolling-origin backtesting of lineage frequency models
#'
#' Evaluates forecast accuracy by repeatedly fitting models on
#' historical data and comparing predictions to held-out observations.
#' This implements the evaluation framework described in Abousamra
#' et al. (2024).
#'
#' @param data An [lfq_data] object.
#' @param engines Character vector of engine names to compare.
#'   Default `"mlr"`.
#' @param origins How to select forecast origins:
#'   * `"weekly"` (default): one origin per unique date, starting
#'     after `min_train` days.
#'   * An integer: use every Nth date as an origin.
#'   * A Date vector: use these specific dates as origins.
#' @param horizons Integer vector of forecast horizons in days.
#'   Default `c(7, 14, 21, 28)`.
#' @param min_train Minimum training window in days. Origins earlier
#'   than `min(date) + min_train` are skipped. Default 42 (6 weeks).
#' @param ... Additional arguments passed to [fit_model()] (e.g.,
#'   `generation_time` for the Piantham engine).
#'
#' @return An `lfq_backtest` object (tibble subclass) with columns:
#' \describe{
#'   \item{origin_date}{Date used as the training cutoff.}
#'   \item{target_date}{Date being predicted.}
#'   \item{horizon}{Forecast horizon in days.}
#'   \item{engine}{Engine name.}
#'   \item{lineage}{Lineage name.}
#'   \item{predicted}{Predicted frequency (median).}
#'   \item{lower}{Lower prediction bound.}
#'   \item{upper}{Upper prediction bound.}
#'   \item{observed}{Observed frequency at target_date.}
#' }
#'
#' @details
#' Implements the rolling-origin evaluation framework described in
#' Abousamra et al. (2024), Section 2.4. At each origin date, the
#' model is fit on data up to that date and forecasts are compared
#' to held-out future observations. This avoids look-ahead bias
#' and provides an honest assessment of real-time forecast accuracy.
#'
#' @references
#' Abousamra E, Figgins M, Bedford T (2024). Fitness models provide
#' accurate short-term forecasts of SARS-CoV-2 variant frequency.
#' \emph{PLoS Computational Biology}, 20(9):e1012443.
#' \doi{10.1371/journal.pcbi.1012443}
#'
#' @seealso [score_forecasts()] to compute accuracy metrics,
#'   [compare_models()] to rank engines.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8),
#'   n_timepoints = 20, seed = 1)
#' bt <- backtest(sim, engines = "mlr",
#'   horizons = c(7, 14), min_train = 42)
#' bt
#' }
#'
#' @export
backtest <- function(data,
                     engines   = "mlr",
                     origins   = "weekly",
                     horizons  = c(7L, 14L, 21L, 28L),
                     min_train = 42L,
                     ...) {

  if (!is_lfq_data(data)) {
    cli::cli_abort("{.arg data} must be an {.cls lfq_data} object.")
  }

  all_dates <- sort(unique(data$.date))
  min_date  <- min(all_dates)
  max_date  <- max(all_dates)
  lineages  <- attr(data, "lineages")

  # --- Determine origin dates ---
  if (identical(origins, "weekly")) {
    eligible <- all_dates[as.numeric(all_dates - min_date) >= min_train]
    if (length(eligible) > 1L) {
      origin_dates <- eligible[-length(eligible)]
    } else {
      cli::cli_abort(
        "Not enough data for backtesting with min_train = {min_train}."
      )
    }
  } else if (is.numeric(origins) && length(origins) == 1L) {
    eligible     <- all_dates[as.numeric(all_dates - min_date) >= min_train]
    origin_dates <- eligible[seq(1L, max(length(eligible) - 1L, 1L),
                                 by = as.integer(origins))]
  } else if (inherits(origins, "Date")) {
    origin_dates <- origins[origins >= min_date + min_train &
                              origins < max_date]
  } else {
    cli::cli_abort(
      "{.arg origins} must be 'weekly', an integer step, or a Date vector."
    )
  }

  if (length(origin_dates) == 0L) {
    cli::cli_abort("No valid origin dates. Try reducing {.arg min_train}.")
  }

  horizons <- as.integer(horizons)

  # --- Main backtest loop ---
  n_total <- length(origin_dates) * length(engines) * length(horizons)
  cli::cli_progress_bar("Backtesting", total = n_total)

  results <- list()

  for (origin in origin_dates) {
    origin <- as.Date(origin, origin = "1970-01-01")

    # Split
    train_data <- data[data$.date <= origin, ]
    train_data <- structure(
      train_data,
      class        = class(data),
      lineages     = attr(data, "lineages"),
      date_range   = c(min(train_data$.date), max(train_data$.date)),
      n_timepoints = length(unique(train_data$.date)),
      has_location = attr(data, "has_location"),
      min_total    = attr(data, "min_total")
    )

    for (eng in engines) {
      # Filter engine-specific args to avoid passing e.g.
      # generation_time to MLR (which doesn't expect it)
      dots <- list(...)
      mlr_args <- c("window", "ci_method", "laplace_smooth")
      piantham_args <- c("generation_time", mlr_args)
      hier_args <- c("shrinkage_method")
      stan_args <- c("chains", "iter_warmup", "iter_sampling")

      allowed <- switch(eng,
        mlr      = mlr_args,
        hier_mlr = c(mlr_args, hier_args),
        piantham = piantham_args,
        fga      = c(stan_args),
        garw     = c(stan_args),
        names(dots)
      )
      engine_dots <- dots[names(dots) %in% allowed]

      fit <- tryCatch(
        do.call(fit_model,
                c(list(data = train_data, engine = eng), engine_dots)),
        error = function(e) NULL
      )

      for (h in horizons) {
        cli::cli_progress_update()
        if (is.null(fit)) next

        fc <- tryCatch(
          forecast(fit, horizon = h, n_sim = 500L),
          error = function(e) NULL
        )
        if (is.null(fc)) next

        target_date <- origin + h

        # Find closest forecast date to target
        fc_forecast <- fc[fc$.type == "forecast", ]
        if (nrow(fc_forecast) == 0L) next

        closest_idx  <- which.min(abs(as.numeric(fc_forecast$.date - target_date)))
        actual_fc_date <- fc_forecast$.date[closest_idx]
        fc_at_target <- fc_forecast[fc_forecast$.date == actual_fc_date, ]

        # Get observed at target (or closest available)
        obs <- data[data$.date == target_date, ]
        if (nrow(obs) == 0L) {
          # Try the closest forecast date instead
          obs <- data[data$.date == actual_fc_date, ]
        }
        if (nrow(obs) == 0L) next

        for (lin in lineages) {
          pred_row <- fc_at_target[fc_at_target$.lineage == lin, ]
          obs_row  <- obs[obs$.lineage == lin, ]

          if (nrow(pred_row) == 1L && nrow(obs_row) == 1L) {
            results <- c(results, list(tibble::tibble(
              origin_date = origin,
              target_date = target_date,
              horizon     = h,
              engine      = eng,
              lineage     = lin,
              predicted   = pred_row$.median,
              lower       = pred_row$.lower,
              upper       = pred_row$.upper,
              observed    = obs_row$.freq
            )))
          }
        }
      }
    }
  }

  cli::cli_progress_done()

  if (length(results) == 0L) {
    cli::cli_warn("Backtesting produced no results.")
    out <- tibble::tibble(
      origin_date = as.Date(character()), target_date = as.Date(character()),
      horizon = integer(), engine = character(), lineage = character(),
      predicted = numeric(), lower = numeric(), upper = numeric(),
      observed = numeric()
    )
  } else {
    out <- dplyr::bind_rows(results)
  }

  structure(
    out,
    class     = c("lfq_backtest", class(tibble::tibble())),
    engines   = engines,
    horizons  = horizons,
    n_origins = length(origin_dates)
  )
}


#' @return The input object, invisibly.
#' @export
print.lfq_backtest <- function(x, ...) {
  n_orig <- attr(x, "n_origins")
  cli::cli_h3("Backtest results")
  cli::cli_text(
    "{nrow(x)} prediction{?s} across {n_orig} origin{?s}"
  )
  cli::cli_text(
    "Engines: {.val {paste(attr(x, 'engines'), collapse = ', ')}}"
  )
  cli::cli_text(
    "Horizons: {paste(attr(x, 'horizons'), collapse = ', ')} days"
  )
  cat("\n")
  NextMethod()
}
