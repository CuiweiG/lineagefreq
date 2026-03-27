#' Piantham engine (internal)
#'
#' Converts MLR growth rates to approximate relative effective
#' reproduction numbers using the Piantham et al. (2021) approach.
#' Does not require case counts; only needs variant frequencies
#' and a generation time.
#'
#' Called by fit_model(engine = "piantham"). Not exported.
#'
#' @noRd
.engine_piantham <- function(data,
                             pivot           = NULL,
                             ci_level        = 0.95,
                             generation_time = NULL,
                             ...) {

  # generation_time is required for this engine
  if (is.null(generation_time) || !is.numeric(generation_time) ||
      length(generation_time) != 1L || generation_time <= 0) {
    cli::cli_abort(
      c("{.arg generation_time} is required for the Piantham engine.",
        "i" = "Pass it via: {.code fit_model(data, engine = 'piantham', generation_time = 5)}")
    )
  }

  # Fit MLR first (all the heavy lifting)
  mlr_result <- .engine_mlr(data, pivot = pivot, ci_level = ci_level, ...)

  # Compute relative Rt for each lineage
  ts         <- mlr_result$time_scale
  tau_scaled <- generation_time / ts   # generation time in time_scale units

  relative_Rt <- exp(mlr_result$growth_rates * tau_scaled)

  # Add Piantham-specific fields
  mlr_result$relative_Rt     <- relative_Rt
  mlr_result$generation_time <- generation_time

  mlr_result
}
