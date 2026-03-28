#' Convert lfq_data to long-format tibble
#'
#' @param x An [lfq_data] object.
#' @param ... Ignored.
#'
#' @return A tibble with all columns.
#'
#' @examples
#' data(sarscov2_us_2022)
#' x <- lfq_data(sarscov2_us_2022, lineage = variant,
#'               date = date, count = count, total = total)
#' as.data.frame(x)
#'
#' @export
as.data.frame.lfq_data <- function(x, ...) {
  tibble::as_tibble(x)
}

#' Convert lfq_fit results to a summary tibble
#'
#' Returns a one-row-per-lineage summary with growth rates,
#' fitted frequencies at first and last time points, and
#' growth advantage in multiple scales.
#'
#' @param fit An `lfq_fit` object.
#' @param generation_time Generation time for Rt calculation.
#'
#' @return A tibble.
#'
#' @examples
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
#' fit <- fit_model(sim)
#' lfq_summary(fit, generation_time = 5)
#'
#' @export
lfq_summary <- function(fit, generation_time = NULL) {
  if (!inherits(fit, "lfq_fit"))
    cli::cli_abort("{.arg fit} must be an {.cls lfq_fit} object.")

  gr <- growth_advantage(fit, type = "growth_rate")

  out <- tibble::tibble(
    lineage     = gr$lineage,
    growth_rate = gr$estimate,
    gr_lower    = gr$lower,
    gr_upper    = gr$upper,
    is_pivot    = gr$lineage == fit$pivot
  )

  if (!is.null(generation_time)) {
    rt <- growth_advantage(fit, type = "relative_Rt",
                           generation_time = generation_time)
    out$relative_Rt <- rt$estimate
    out$Rt_lower    <- rt$lower
    out$Rt_upper    <- rt$upper
  }

  out
}
