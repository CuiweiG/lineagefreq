#' Plot lineage frequency model results
#'
#' @param object An `lfq_fit` object.
#' @param type Plot type: `"frequency"` (default), `"advantage"`,
#'   `"trajectory"`, or `"residuals"`.
#' @param generation_time Required when `type = "advantage"`.
#' @param ... Ignored.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
#' fit <- fit_model(sim)
#' autoplot(fit)
#' autoplot(fit, type = "advantage", generation_time = 5)
#' }
#'
#' @export
autoplot.lfq_fit <- function(object,
                             type = c("frequency", "advantage",
                                      "trajectory", "residuals"),
                             generation_time = NULL, ...) {
  type <- match.arg(type)
  switch(type,
    frequency  = .plot_frequency(object),
    advantage  = .plot_advantage(object, generation_time),
    trajectory = .plot_trajectory(object),
    residuals  = .plot_residuals(object)
  )
}


#' Plot a lineage frequency forecast
#'
#' @param object An `lfq_forecast` object.
#' @param ... Ignored.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
#' fit <- fit_model(sim)
#' fc <- forecast(fit, horizon = 14)
#' autoplot(fc)
#' }
#'
#' @export
autoplot.lfq_forecast <- function(object, ...) {
  fc_start <- as.numeric(as.Date(attr(object, "forecast_start")))
  lins     <- attr(object, "lineages")

  ggplot2::ggplot(object, ggplot2::aes(
    x = .data$.date, y = .data$.median,
    colour = .data$.lineage, fill = .data$.lineage
  )) +
    ggplot2::geom_line(
      data = function(d) d[d$.type == "fitted", ],
      linewidth = 0.8
    ) +
    ggplot2::geom_ribbon(
      data = function(d) d[d$.type == "forecast", ],
      ggplot2::aes(ymin = .data$.lower, ymax = .data$.upper),
      alpha = 0.2, colour = NA
    ) +
    ggplot2::geom_line(
      data = function(d) d[d$.type == "forecast", ],
      linetype = "dashed", linewidth = 0.6
    ) +
    ggplot2::geom_vline(
      xintercept = fc_start, linetype = "dotted", colour = "gray40"
    ) +
    ggplot2::scale_x_date(date_labels = "%Y-%m-%d") +
    ggplot2::scale_y_continuous(
      labels = .pct_label, limits = c(0, 1), expand = c(0, 0)
    ) +
    .lfq_scale_colour(lins) +
    .lfq_scale_fill(lins) +
    ggplot2::labs(x = NULL, y = "Frequency",
                  colour = "Lineage", fill = "Lineage") +
    .lfq_theme()
}


#  Internal helpers 

#' @noRd
.plot_frequency <- function(fit) {
  df   <- fit$residuals
  lins <- fit$lineages

  ggplot2::ggplot(df, ggplot2::aes(x = .data$.date)) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$.observed, colour = .data$.lineage),
      size = 1.8, alpha = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$.fitted_freq, colour = .data$.lineage),
      linewidth = 1.0
    ) +
    ggplot2::scale_x_date(date_labels = "%Y-%m-%d") +
    ggplot2::scale_y_continuous(
      labels = .pct_label, limits = c(0, 1), expand = c(0.01, 0)
    ) +
    .lfq_scale_colour(lins) +
    ggplot2::labs(x = NULL, y = "Frequency", colour = "Lineage") +
    .lfq_theme()
}


#' @noRd
.plot_advantage <- function(fit, generation_time) {
  if (is.null(generation_time)) {
    cli::cli_abort(
      "{.arg generation_time} is required for type {.val advantage}."
    )
  }

  ga   <- growth_advantage(fit, type = "relative_Rt",
                           generation_time = generation_time)
  ga   <- ga[ga$lineage != fit$pivot, ]
  lins <- ga$lineage

  ggplot2::ggplot(ga, ggplot2::aes(
    x    = .data$estimate,
    y    = stats::reorder(.data$lineage, .data$estimate),
    xmin = .data$lower, xmax = .data$upper,
    colour = .data$lineage
  )) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                        colour = "gray50") +
    ggplot2::geom_pointrange(linewidth = 0.6, size = 0.8) +
    .lfq_scale_colour(lins) +
    ggplot2::labs(
      x       = paste0("Relative Rt (gen. time = ",
                       generation_time, " days)"),
      y       = NULL,
      caption = paste("Pivot:", fit$pivot)
    ) +
    .lfq_theme() +
    ggplot2::theme(legend.position = "none")
}


#' @noRd
.plot_trajectory <- function(fit) {
  df   <- fit$residuals
  lins <- fit$lineages

  ggplot2::ggplot(df, ggplot2::aes(x = .data$.date)) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$.observed, colour = .data$.lineage),
      alpha = 0.5, size = 1.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$.fitted_freq, colour = .data$.lineage),
      linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~ .data$.lineage, scales = "free_y") +
    ggplot2::scale_x_date(date_labels = "%Y-%m") +
    ggplot2::scale_y_continuous(labels = .pct_label) +
    .lfq_scale_colour(lins) +
    ggplot2::labs(x = NULL, y = "Frequency") +
    .lfq_theme(base_size = 11) +
    ggplot2::theme(legend.position = "none")
}


#' @noRd
.plot_residuals <- function(fit) {
  df   <- fit$residuals
  df   <- df[is.finite(df$.pearson_resid), ]
  lins <- fit$lineages

  ggplot2::ggplot(df, ggplot2::aes(
    x = .data$.date, y = .data$.pearson_resid,
    colour = .data$.lineage
  )) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_x_date(date_labels = "%Y-%m-%d") +
    .lfq_scale_colour(lins) +
    ggplot2::labs(x = NULL, y = "Pearson residual",
                  colour = "Lineage") +
    .lfq_theme()
}


#' Percent formatter (no scales dependency)
#' @noRd
.pct_label <- function(x) paste0(round(x * 100), "%")
