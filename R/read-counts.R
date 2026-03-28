#' Read lineage count data from a CSV file
#'
#' A convenience wrapper for reading surveillance count data from
#' CSV files into [lfq_data] format. Expects a file with at least
#' columns for date, lineage, and count.
#'
#' @param file Path to CSV file.
#' @param date Name of the date column. Default `"date"`.
#' @param lineage Name of the lineage column. Default `"lineage"`.
#' @param count Name of the count column. Default `"count"`.
#' @param ... Additional arguments passed to [lfq_data()].
#'
#' @return An [lfq_data] object.
#'
#' @examples
#' # Read the bundled example CSV
#' f <- system.file("extdata", "example_counts.csv",
#'                  package = "lineagefreq")
#' x <- read_lineage_counts(f)
#' x
#'
#' @export
read_lineage_counts <- function(file,
                                date    = "date",
                                lineage = "lineage",
                                count   = "count",
                                ...) {
  df <- utils::read.csv(file, stringsAsFactors = FALSE)
  if (!date %in% names(df))
    cli::cli_abort("Column {.field {date}} not found in {.arg file}.")
  if (!lineage %in% names(df))
    cli::cli_abort("Column {.field {lineage}} not found in {.arg file}.")
  if (!count %in% names(df))
    cli::cli_abort("Column {.field {count}} not found in {.arg file}.")

  df[[date]] <- as.Date(df[[date]])
  lfq_data(df,
           lineage = !!rlang::sym(lineage),
           date    = !!rlang::sym(date),
           count   = !!rlang::sym(count),
           ...)
}
