#' Create a lineage frequency data object
#'
#' Validates, structures, and annotates lineage count data for
#' downstream modeling and analysis. This is the entry point for
#' all lineagefreq workflows.
#'
#' @param data A data frame containing at minimum columns for lineage
#'   identity, date, and count.
#' @param lineage <[`tidy-select`][dplyr::dplyr_tidy_select]> Column
#'   containing lineage/variant identifiers (character or factor).
#' @param date <[`tidy-select`][dplyr::dplyr_tidy_select]> Column
#'   containing collection dates (`Date` class or parseable character).
#' @param count <[`tidy-select`][dplyr::dplyr_tidy_select]> Column
#'   containing sequence counts (non-negative integers).
#' @param total <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional
#'   column of total sequences per date-location. If `NULL`, computed
#'   as the sum of `count` per group.
#' @param location <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional
#'   column for geographic stratification.
#' @param min_total Minimum total count per time point. Time points
#'   below this are flagged as unreliable. Default 10.
#'
#' @return An `lfq_data` object (a tibble subclass) with standardized
#'   columns:
#'   \describe{
#'     \item{`.lineage`}{Lineage identifier (character).}
#'     \item{`.date`}{Collection date (Date).}
#'     \item{`.count`}{Sequence count (integer).}
#'     \item{`.total`}{Total sequences at this time point (integer).}
#'     \item{`.freq`}{Observed frequency (numeric).}
#'     \item{`.reliable`}{Logical; `TRUE` if `.total >= min_total`.}
#'     \item{`.location`}{Location, if provided (character).}
#'   }
#'   All original columns are preserved.
#'
#' @details
#' Performs the following validation and processing:
#' \enumerate{
#'   \item Checks that all required columns exist and have correct types.
#'   \item Coerces character dates to Date via ISO 8601 parsing.
#'   \item Ensures counts are non-negative integers.
#'   \item Replaces NA counts with 0 (with warning).
#'   \item Aggregates duplicate lineage-date rows by summing (with warning).
#'   \item Computes per-time-point totals and frequencies.
#'   \item Flags time points below `min_total` as unreliable.
#'   \item Sorts by date ascending, then lineage alphabetically.
#' }
#'
#' @examples
#' d <- data.frame(
#'   date = rep(seq(as.Date("2024-01-01"), by = "week",
#'                  length.out = 8), each = 3),
#'   lineage = rep(c("JN.1", "KP.3", "Other"), 8),
#'   n = c(5, 2, 93, 12, 5, 83, 28, 11, 61, 50, 20, 30,
#'         68, 18, 14, 80, 12, 8, 88, 8, 4, 92, 5, 3)
#' )
#' x <- lfq_data(d, lineage = lineage, date = date, count = n)
#' x
#'
#' @export
lfq_data <- function(data, lineage, date, count,
                     total = NULL, location = NULL,
                     min_total = 10L) {

  # --- Capture column names via tidy evaluation ---
  lin_col  <- rlang::as_name(rlang::enquo(lineage))
  date_col <- rlang::as_name(rlang::enquo(date))
  cnt_col  <- rlang::as_name(rlang::enquo(count))

  tot_col <- NULL
  if (!rlang::quo_is_null(rlang::enquo(total))) {
    tot_col <- rlang::as_name(rlang::enquo(total))
  }

  loc_col <- NULL
  if (!rlang::quo_is_null(rlang::enquo(location))) {
    loc_col <- rlang::as_name(rlang::enquo(location))
  }

  # --- Logic 1: Check required columns exist ---
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }
  assert_col(data, lin_col,  "data")
  assert_col(data, date_col, "data")
  assert_col(data, cnt_col,  "data")
  if (!is.null(tot_col)) assert_col(data, tot_col, "data")
  if (!is.null(loc_col)) assert_col(data, loc_col, "data")

  df <- tibble::as_tibble(data)

  # --- Logic 2: Coerce date ---
  if (is.character(df[[date_col]])) {
    parsed <- tryCatch(as.Date(df[[date_col]]), error = function(e) NULL)
    if (is.null(parsed) || all(is.na(parsed))) {
      cli::cli_abort(
        "Column {.field {date_col}} could not be parsed as Date."
      )
    }
    df[[date_col]] <- parsed
  }
  assert_date(df[[date_col]], date_col)

  # --- Logic 3 & 4: Validate and coerce count ---
  cnt <- df[[cnt_col]]
  if (!is.numeric(cnt)) {
    cli::cli_abort(
      "Column {.field {cnt_col}} must be numeric. Got {.cls {class(cnt)}}."
    )
  }
  if (any(is.na(cnt))) {
    n_na <- sum(is.na(cnt))
    cli::cli_warn(
      "Replacing {n_na} NA value{?s} in {.field {cnt_col}} with 0."
    )
    cnt[is.na(cnt)] <- 0
  }
  if (any(cnt < 0)) {
    cli::cli_abort(
      "Column {.field {cnt_col}} contains negative values."
    )
  }
  df[[cnt_col]] <- as.integer(round(cnt))

  # --- Coerce lineage to character ---
  df[[lin_col]] <- as.character(df[[lin_col]])

  # --- Logic 5: Aggregate duplicates ---
  grp <- c(date_col, lin_col)
  if (!is.null(loc_col)) grp <- c(grp, loc_col)

  dup_check <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::ungroup()

  if (nrow(dup_check) > 0) {
    n_dup <- nrow(dup_check)
    cli::cli_warn(
      c("!" = "{n_dup} duplicate row{?s} found.",
        "i" = "Aggregating by summing counts.")
    )
    df <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(cnt_col), sum),
        .groups = "drop"
      )
  }

  # --- Logic 6: Compute totals ---
  tgrp <- date_col
  if (!is.null(loc_col)) tgrp <- c(tgrp, loc_col)

  if (is.null(tot_col)) {
    df <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(tgrp))) |>
      dplyr::mutate(.total = sum(.data[[cnt_col]])) |>
      dplyr::ungroup()
  } else {
    df$.total <- as.integer(df[[tot_col]])
  }

  # --- Compute frequency ---
  df$.freq <- ifelse(
    df$.total > 0,
    df[[cnt_col]] / df$.total,
    NA_real_
  )

  # --- Rename to standardized names ---
  df <- dplyr::rename(
    df,
    .lineage = dplyr::all_of(lin_col),
    .date    = dplyr::all_of(date_col),
    .count   = dplyr::all_of(cnt_col)
  )
  if (!is.null(loc_col)) {
    df <- dplyr::rename(df, .location = dplyr::all_of(loc_col))
  }

  # --- Logic 7: Reliability flag ---
  df$.reliable <- df$.total >= min_total

  n_unrel <- sum(!df$.reliable) / max(length(unique(df$.lineage)), 1L)
  if (n_unrel > 0) {
    cli::cli_inform(
      "{.val {as.integer(n_unrel)}} time point{?s} have fewer than {min_total} total sequences."
    )
  }

  # --- Logic 8: Sort ---
  sort_cols <- c(".date", ".lineage")
  if (!is.null(loc_col)) sort_cols <- c(".location", sort_cols)
  df <- dplyr::arrange(df, dplyr::across(dplyr::all_of(sort_cols)))

  # --- Construct S3 class ---
  structure(
    df,
    class       = c("lfq_data", class(tibble::tibble())),
    lineages    = sort(unique(df$.lineage)),
    date_range  = range(df$.date, na.rm = TRUE),
    n_timepoints = length(unique(df$.date)),
    has_location = !is.null(loc_col),
    min_total   = min_total
  )
}


#' Test if an object is an lfq_data object
#'
#' @param x Object to test.
#'
#' @return Logical scalar.
#'
#' @examples
#' d <- data.frame(date = Sys.Date(), lineage = "A", count = 10)
#' x <- lfq_data(d, lineage = lineage, date = date, count = count)
#' is_lfq_data(x)
#' is_lfq_data(d)
#'
#' @export
is_lfq_data <- function(x) {
  inherits(x, "lfq_data")
}


#' @export
print.lfq_data <- function(x, ...) {
  n_lin <- length(attr(x, "lineages"))
  n_tp  <- attr(x, "n_timepoints")
  dr    <- attr(x, "date_range")

  cli::cli_h3("Lineage frequency data")
  cli::cli_text(
    "{.val {n_lin}} lineage{?s}, {.val {n_tp}} time point{?s}"
  )
  cli::cli_text("Date range: {dr[1]} to {dr[2]}")

  lin_display <- utils::head(attr(x, "lineages"), 5)
  cli::cli_text(
    "Lineages: {.val {paste(lin_display, collapse = ', ')}}{if (n_lin > 5) ', ...' else ''}"
  )

  n_unrel <- sum(!x$.reliable) / max(n_lin, 1L)
  if (n_unrel > 0) {
    cli::cli_text(
      "{.val {as.integer(n_unrel)}} time point{?s} flagged as low-count"
    )
  }

  cat("\n")
  NextMethod()
}
