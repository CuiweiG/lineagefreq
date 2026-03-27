#' Summarize emerging lineages
#'
#' Tests whether each lineage's frequency is significantly increasing
#' over time using a binomial GLM. Useful for early warning of
#' lineages that may warrant enhanced surveillance.
#'
#' @param data An [lfq_data] object.
#' @param threshold Minimum current frequency to test. Default 0.01.
#' @param min_obs Minimum time points observed. Default 3.
#' @param p_adjust P-value adjustment method. Default `"holm"`.
#'
#' @return A tibble with columns: `lineage`, `first_seen`, `last_seen`,
#'   `n_timepoints`, `current_freq`, `growth_rate`, `p_value`,
#'   `p_adjusted`, `significant`, `direction`.
#'
#' @examples
#' sim <- simulate_dynamics(
#'   n_lineages = 4,
#'   advantages = c(emerging = 1.5, stable = 1.0, declining = 0.7),
#'   n_timepoints = 12, seed = 42)
#' summarize_emerging(sim)
#'
#' @export
summarize_emerging <- function(data,
                               threshold = 0.01,
                               min_obs   = 3L,
                               p_adjust  = "holm") {

  if (!is_lfq_data(data)) {
    cli::cli_abort("{.arg data} must be an {.cls lfq_data} object.")
  }

  lineages <- attr(data, "lineages")
  results  <- list()

  for (lin in lineages) {
    lin_data <- data[data$.lineage == lin & data$.reliable, ]
    n_tp     <- nrow(lin_data)
    if (n_tp < min_obs) next

    current_freq <- utils::tail(lin_data$.freq, 1L)
    if (is.na(current_freq) || current_freq < threshold) next

    t_num <- as.numeric(lin_data$.date - min(lin_data$.date))

    fit_full <- tryCatch(
      stats::glm(
        cbind(.count, .total - .count) ~ t_num,
        data   = lin_data,
        family = stats::binomial()
      ),
      error = function(e) NULL
    )
    if (is.null(fit_full)) next

    fit_null <- tryCatch(
      stats::glm(
        cbind(.count, .total - .count) ~ 1,
        data   = lin_data,
        family = stats::binomial()
      ),
      error = function(e) NULL
    )
    if (is.null(fit_null)) next

    lr_stat   <- -2 * (as.numeric(stats::logLik(fit_null)) -
                          as.numeric(stats::logLik(fit_full)))
    p_val     <- stats::pchisq(lr_stat, df = 1L, lower.tail = FALSE)
    beta_time <- unname(stats::coef(fit_full)["t_num"])

    results <- c(results, list(tibble::tibble(
      lineage      = lin,
      first_seen   = min(lin_data$.date),
      last_seen    = max(lin_data$.date),
      n_timepoints = n_tp,
      current_freq = current_freq,
      growth_rate  = beta_time,
      p_value      = p_val
    )))
  }

  if (length(results) == 0L) {
    return(tibble::tibble(
      lineage = character(), first_seen = as.Date(character()),
      last_seen = as.Date(character()), n_timepoints = integer(),
      current_freq = numeric(), growth_rate = numeric(),
      p_value = numeric(), p_adjusted = numeric(),
      significant = logical(), direction = character()
    ))
  }

  out <- dplyr::bind_rows(results)
  out$p_adjusted  <- stats::p.adjust(out$p_value, method = p_adjust)
  out$significant <- out$p_adjusted < 0.05
  out$direction   <- dplyr::case_when(
    out$growth_rate >  0.001 & out$significant ~ "growing",
    out$growth_rate < -0.001 & out$significant ~ "declining",
    TRUE ~ "stable"
  )
  dplyr::arrange(out, .data$p_adjusted)
}
