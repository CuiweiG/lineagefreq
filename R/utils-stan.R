#' Get path to a bundled Stan model file
#' @noRd
.stan_model_path <- function(model_name) {
  system.file("stan", paste0(model_name, ".stan"),
              package = "lineagefreq", mustWork = TRUE)
}

#' Compile a Stan model (with caching)
#' @noRd
.get_stan_model <- function(model_name) {
  rlang::check_installed("cmdstanr",
    reason = "to fit Bayesian models (FGA/GARW).",
    action = function(pkg, ...) {
      message(
        'Install with: install.packages("cmdstanr", ',
        'repos = c("https://mc-stan.org/r-packages/", getOption("repos")))'
      )
    }
  )

  path <- .stan_model_path(model_name)
  cmdstanr::cmdstan_model(path, quiet = TRUE)
}

#' Prepare Stan data from lfq_data
#' @noRd
.prepare_stan_data <- function(data, pivot) {

  lineages <- attr(data, "lineages")
  dates    <- sort(unique(data$.date))
  n_tp     <- length(dates)
  n_lin    <- length(lineages)

  # Select pivot
  if (is.null(pivot)) {
    first_date <- dates[1L]
    early      <- data[data$.date == first_date, ]
    totals_by_lin <- stats::setNames(early$.count, early$.lineage)
    pivot <- names(which.max(totals_by_lin))
  }
  pivot_idx <- match(pivot, lineages)

  # Reshape to count matrix
  wide <- data[, c(".date", ".lineage", ".count")]
  wide <- tidyr::pivot_wider(
    wide,
    names_from  = ".lineage",
    values_from = ".count",
    values_fill = 0L
  )
  wide <- dplyr::arrange(wide, .data$.date)
  Y    <- as.matrix(wide[, lineages, drop = FALSE])

  list(
    T            = n_tp,
    V            = n_lin,
    Y            = Y,
    pivot        = pivot_idx,
    lineages     = lineages,
    pivot_name   = pivot,
    dates        = dates,
    n_timepoints = n_tp
  )
}
