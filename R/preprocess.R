#' Collapse rare lineages into an aggregate group
#'
#' Merges lineages that never exceed a frequency or count threshold
#' into a single group (default "Other"). Useful for reducing noise
#' from dozens of low-frequency lineages.
#'
#' @param x An [lfq_data] object.
#' @param min_freq Minimum peak frequency a lineage must reach at
#'   any time point to be kept. Default 0.01 (1%).
#' @param min_count Minimum total count across all time points to
#'   be kept. Default 10.
#' @param other_label Label for the collapsed group. Default "Other".
#' @param mapping Optional named character vector for custom grouping.
#'   Names = original lineage names, values = group names. If provided,
#'   `min_freq` and `min_count` are ignored.
#'
#' @return An [lfq_data] object with rare lineages merged.
#'
#' @examples
#' sim <- simulate_dynamics(n_lineages = 6,
#'   advantages = c(A = 1.3, B = 1.1, C = 0.95, D = 0.5, E = 0.3),
#'   n_timepoints = 12, seed = 1)
#' collapsed <- collapse_lineages(sim, min_freq = 0.05)
#' attr(collapsed, "lineages")
#'
#' @export
collapse_lineages <- function(x,
                              min_freq    = 0.01,
                              min_count   = 10L,
                              other_label = "Other",
                              mapping     = NULL) {

  if (!is_lfq_data(x)) {
    cli::cli_abort("{.arg x} must be an {.cls lfq_data} object.")
  }

  df <- tibble::as_tibble(x)

  if (!is.null(mapping)) {
    df$.lineage <- ifelse(
      df$.lineage %in% names(mapping),
      mapping[df$.lineage],
      df$.lineage
    )
  } else {
    lin_stats <- df |>
      dplyr::group_by(.data$.lineage) |>
      dplyr::summarise(
        peak_freq   = max(.data$.freq, na.rm = TRUE),
        total_count = sum(.data$.count, na.rm = TRUE),
        .groups     = "drop"
      )

    rare <- lin_stats$.lineage[
      lin_stats$peak_freq < min_freq | lin_stats$total_count < min_count
    ]

    if (length(rare) > 0L) {
      cli::cli_inform(
        "Collapsing {length(rare)} rare lineage{?s} into {.val {other_label}}."
      )
      df$.lineage[df$.lineage %in% rare] <- other_label
    }
  }

  # Re-aggregate
  grp <- ".date"
  if (".location" %in% names(df)) grp <- c(grp, ".location")
  grp <- c(grp, ".lineage")

  df <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
    dplyr::summarise(.count = sum(.data$.count), .groups = "drop")

  total_grp <- grp[grp != ".lineage"]
  df <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(total_grp))) |>
    dplyr::mutate(
      .total = sum(.data$.count),
      .freq  = ifelse(.data$.total > 0, .data$.count / .data$.total,
                      NA_real_)
    ) |>
    dplyr::ungroup()

  df$.reliable <- df$.total >= attr(x, "min_total")
  df <- dplyr::arrange(df, dplyr::across(dplyr::all_of(grp)))

  structure(
    df,
    class        = c("lfq_data", class(tibble::tibble())),
    lineages     = sort(unique(df$.lineage)),
    date_range   = attr(x, "date_range"),
    n_timepoints = attr(x, "n_timepoints"),
    has_location = attr(x, "has_location"),
    min_total    = attr(x, "min_total")
  )
}


#' Filter sparse time points and lineages
#'
#' Removes time points with very low total counts and lineages
#' observed at too few time points.
#'
#' @param x An [lfq_data] object.
#' @param min_total Minimum total sequences per time point. Default 10.
#' @param min_timepoints Minimum number of time points a lineage
#'   must appear to be retained. Default 3.
#'
#' @return An [lfq_data] object with sparse entries removed.
#'
#' @examples
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.2, "B" = 0.8), seed = 1)
#' filtered <- filter_sparse(sim, min_total = 100)
#'
#' @export
filter_sparse <- function(x, min_total = 10L, min_timepoints = 3L) {

  if (!is_lfq_data(x)) {
    cli::cli_abort("{.arg x} must be an {.cls lfq_data} object.")
  }

  df <- tibble::as_tibble(x)

  # Filter by total
  df <- df[df$.total >= min_total, ]

  # Filter lineages by time points present
  lin_tp  <- df |>
    dplyr::group_by(.data$.lineage) |>
    dplyr::summarise(n_tp = dplyr::n_distinct(.data$.date),
                     .groups = "drop")
  keep_lin <- lin_tp$.lineage[lin_tp$n_tp >= min_timepoints]
  df <- df[df$.lineage %in% keep_lin, ]

  if (nrow(df) == 0L) {
    cli::cli_warn(
      "All data removed by filtering. Try less strict thresholds."
    )
  }

  structure(
    df,
    class        = c("lfq_data", class(tibble::tibble())),
    lineages     = sort(unique(df$.lineage)),
    date_range   = if (nrow(df) > 0L) range(df$.date) else as.Date(c(NA, NA)),
    n_timepoints = length(unique(df$.date)),
    has_location = attr(x, "has_location"),
    min_total    = min_total
  )
}
