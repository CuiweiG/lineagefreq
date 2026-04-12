#' Construct a population immunity landscape
#'
#' Assembles time-varying population immunity estimates against
#' each circulating lineage from seroprevalence surveys, vaccination
#' records, or infection history data. The resulting object serves
#' as input to \code{\link{fitness_decomposition}} for disentangling
#' intrinsic transmissibility from immune escape.
#'
#' @param data A data frame with columns for date, lineage, and
#'   immunity level.
#' @param date Column name (unquoted or string) containing dates.
#' @param lineage Column name containing lineage/variant identifiers.
#' @param immunity Column name containing population-level immunity
#'   estimates (proportion, 0--1 scale) against each lineage.
#' @param type Character vector specifying immunity source(s):
#'   \code{"infection"}, \code{"vaccine"}, \code{"hybrid"}, or
#'   \code{"combined"} (default). If the data contain a \code{type}
#'   column, per-source estimates are preserved; otherwise all
#'   estimates are treated as combined.
#' @param cross_immunity Optional numeric matrix of cross-immunity
#'   between lineages. Rows and columns correspond to lineages;
#'   entry \eqn{(i,j)} is the degree to which immunity against
#'   lineage \eqn{i} protects against lineage \eqn{j} (0 to 1).
#'   Default \code{NULL} assumes no cross-protection data.
#'
#' @return An \code{immune_landscape} object (S3 class) with
#'   components:
#'   \describe{
#'     \item{estimates}{Tibble with columns \code{date}, \code{lineage},
#'       \code{immunity}, and optionally \code{type}.}
#'     \item{lineages}{Character vector of lineage names.}
#'     \item{date_range}{Two-element Date vector.}
#'     \item{cross_immunity}{Matrix or NULL.}
#'   }
#'
#' @details
#' Immunity estimates may come from multiple sources. Seroprevalence
#' surveys provide direct measurements but are infrequent and
#' geographically sparse. Vaccination coverage data are more widely
#' available but do not capture waning or variant-specific escape.
#' Model-based reconstructions (e.g., from case and death data) can
#' fill gaps but introduce model dependence.
#'
#' \code{immune_landscape()} accepts any of these as input. The
#' critical requirement is that immunity is expressed on a 0--1
#' scale representing the proportion of the population with
#' neutralising protection against each lineage at each time point.
#'
#' @seealso \code{\link{fitness_decomposition}} for downstream
#'   analysis.
#'
#' @examples
#' # Simulated immunity data
#' imm_data <- data.frame(
#'   date = rep(seq(as.Date("2024-01-01"), by = "week",
#'                  length.out = 10), each = 3),
#'   lineage = rep(c("BA.5", "XBB.1.5", "JN.1"), 10),
#'   immunity = c(
#'     rep(c(0.6, 0.4, 0.1), 5),
#'     rep(c(0.55, 0.45, 0.25), 5))
#' )
#' il <- immune_landscape(imm_data, date = date,
#'   lineage = lineage, immunity = immunity)
#' il
#'
#' @export
immune_landscape <- function(data,
                             date,
                             lineage,
                             immunity,
                             type = "combined",
                             cross_immunity = NULL) {

  date_col <- as.character(substitute(date))
  lin_col  <- as.character(substitute(lineage))
  imm_col  <- as.character(substitute(immunity))

  # Handle string or symbol column references
  if (!date_col %in% names(data)) {
    date_col <- deparse(substitute(date))
  }
  if (!lin_col %in% names(data)) {
    lin_col <- deparse(substitute(lineage))
  }
  if (!imm_col %in% names(data)) {
    imm_col <- deparse(substitute(immunity))
  }

  assert_col(data, date_col, "date")
  assert_col(data, lin_col, "lineage")
  assert_col(data, imm_col, "immunity")

  estimates <- tibble::tibble(
    date     = as.Date(data[[date_col]]),
    lineage  = as.character(data[[lin_col]]),
    immunity = as.numeric(data[[imm_col]])
  )

  if ("type" %in% names(data)) {
    estimates$type <- as.character(data[["type"]])
  } else {
    estimates$type <- type[1]
  }

  # Validate immunity range
  if (any(estimates$immunity < 0 | estimates$immunity > 1,
          na.rm = TRUE)) {
    cli::cli_abort("Immunity values must be between 0 and 1.")
  }

  lineages <- sort(unique(estimates$lineage))

  # Validate cross-immunity matrix if provided
  if (!is.null(cross_immunity)) {
    if (!is.matrix(cross_immunity)) {
      cli::cli_abort("{.arg cross_immunity} must be a numeric matrix.")
    }
    if (nrow(cross_immunity) != length(lineages) ||
        ncol(cross_immunity) != length(lineages)) {
      cli::cli_abort(
        "{.arg cross_immunity} dimensions must match number of lineages ({length(lineages)})."
      )
    }
    rownames(cross_immunity) <- lineages
    colnames(cross_immunity) <- lineages
  }

  structure(
    list(
      estimates      = estimates,
      lineages       = lineages,
      date_range     = range(estimates$date),
      cross_immunity = cross_immunity
    ),
    class = "immune_landscape"
  )
}


#' @export
print.immune_landscape <- function(x, ...) {
  cli::cli_h3("Population immunity landscape")
  cli::cli_text("{length(x$lineages)} lineage{?s}: {.val {paste(x$lineages, collapse = ', ')}}")
  cli::cli_text("Date range: {x$date_range[1]} to {x$date_range[2]}")
  types <- unique(x$estimates$type)
  cli::cli_text("Immunity types: {.val {paste(types, collapse = ', ')}}")
  if (!is.null(x$cross_immunity)) {
    cli::cli_text("Cross-immunity matrix: provided")
  }
  invisible(x)
}


#' Plot population immunity landscape
#'
#' @param x An \code{immune_landscape} object.
#' @param ... Unused.
#' @return A ggplot object.
#' @export
plot.immune_landscape <- function(x, ...) {
  ggplot2::ggplot(x$estimates,
    ggplot2::aes(x = .data$date, y = .data$immunity,
                 colour = .data$lineage)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Date", y = "Population immunity",
      colour = "Lineage",
      title = "Population immunity landscape"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}
