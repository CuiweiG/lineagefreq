#' Coerce to lfq_data
#'
#' Generic function to convert various data formats into [lfq_data]
#' objects. Methods can be defined for specific input classes to enable
#' seamless interoperability with other genomic surveillance packages.
#'
#' @param x An object to coerce.
#' @param ... Additional arguments passed to methods. For the
#'   `data.frame` method, these are passed to [lfq_data()].
#'
#' @return An [lfq_data] object.
#'
#' @seealso [lfq_data()] for the primary constructor.
#'
#' @examples
#' df <- data.frame(
#'   date    = rep(as.Date("2024-01-01") + c(0, 7), each = 2),
#'   lineage = rep(c("A", "B"), 2),
#'   count   = c(80, 20, 60, 40)
#' )
#' x <- as_lfq_data(df, lineage = lineage, date = date, count = count)
#' x
#'
#' @export
as_lfq_data <- function(x, ...) {
  UseMethod("as_lfq_data")
}

#' @rdname as_lfq_data
#' @return An [lfq_data] object.
#' @export
as_lfq_data.lfq_data <- function(x, ...) {
  x
}

#' @rdname as_lfq_data
#' @return An [lfq_data] object.
#' @export
as_lfq_data.data.frame <- function(x, ...) {
  lfq_data(x, ...)
}
