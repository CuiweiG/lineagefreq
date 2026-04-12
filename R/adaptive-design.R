#' Adaptive sequencing allocation via Thompson sampling
#'
#' Reallocates sequencing resources across regions in real time,
#' directing effort toward strata where uncertainty reduction has
#' the highest decision value. Unlike static Neyman allocation
#' (which optimises for a single time point), this function adapts
#' to evolving variant dynamics across multiple surveillance rounds.
#'
#' @param data An \code{lfq_data} object with a location column.
#' @param capacity Total sequencing capacity per round (integer).
#' @param n_rounds Number of allocation rounds to simulate. Default
#'   \code{NULL} uses the number of time points in the data.
#' @param strategy Allocation strategy: \code{"thompson"} (default)
#'   for Thompson sampling, or \code{"ucb"} for Upper Confidence
#'   Bound.
#' @param target_lineage Character; lineage to optimise detection
#'   for. Default \code{NULL} optimises for overall frequency
#'   estimation.
#' @param exploration Exploration parameter for UCB. Default 2.0.
#'   Larger values explore more; smaller values exploit current
#'   best regions.
#' @param seed Random seed. Default \code{NULL}.
#'
#' @return An \code{adaptive_allocation} S3 class with components:
#'   \describe{
#'     \item{allocations}{Tibble with \code{round}, \code{region},
#'       \code{n_allocated}, \code{uncertainty}, \code{frequency}.}
#'     \item{summary}{Tibble with per-region totals and mean
#'       allocation.}
#'     \item{strategy}{Character; strategy used.}
#'     \item{capacity}{Integer; per-round capacity.}
#'   }
#'
#' @details
#' \strong{Thompson sampling} draws from the posterior distribution
#' of frequency estimates and allocates proportional to the
#' posterior variance. Regions with high uncertainty receive more
#' sequences, but the stochastic draws naturally balance
#' exploration (sampling uncertain regions) and exploitation
#' (sampling where variants are most prevalent).
#'
#' \strong{UCB} allocates proportional to the upper confidence bound
#' of the estimation error: \eqn{\text{score}_r = \hat{\sigma}_r +
#' c \sqrt{2 \log(t) / n_r}} where \eqn{\hat{\sigma}_r} is the
#' current estimation uncertainty in region \eqn{r}, \eqn{n_r} is
#' the cumulative allocation, \eqn{t} is the round, and \eqn{c}
#' is the exploration parameter.
#'
#' @seealso \code{\link{surveillance_value}} for EVOI analysis.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.3, "B" = 0.9),
#'   n_timepoints = 12, seed = 1)
#' ad <- adaptive_design(sim, capacity = 200, n_rounds = 8)
#' ad
#' }
#'
#' @export
adaptive_design <- function(data,
                            capacity,
                            n_rounds       = NULL,
                            strategy       = c("thompson", "ucb"),
                            target_lineage = NULL,
                            exploration    = 2.0,
                            seed           = NULL) {

  if (!is_lfq_data(data)) {
    cli::cli_abort("{.arg data} must be an {.cls lfq_data} object.")
  }
  assert_pos_int(capacity, "capacity")
  strategy <- match.arg(strategy)
  if (!is.null(seed)) set.seed(seed)

  dates <- sort(unique(data$.date))
  if (is.null(n_rounds)) n_rounds <- length(dates)
  n_rounds <- min(n_rounds, length(dates))

  lineages <- attr(data, "lineages")

  # Use lineages as "regions" for allocation if no location column
  # In practice, regions would come from a location column
  has_loc <- attr(data, "has_location")
  if (has_loc && ".location" %in% names(data)) {
    regions <- sort(unique(data$.location))
  } else {
    # Treat each time point as an allocation decision
    regions <- paste0("region_", seq_len(min(5L, n_rounds)))
  }

  n_regions <- length(regions)
  cum_alloc <- rep(1L, n_regions)  # initialise with 1 to avoid division by zero
  names(cum_alloc) <- regions

  alloc_rows <- list()

  for (round in seq_len(n_rounds)) {
    # Estimate current uncertainty per region
    # Use variance of observed frequencies as proxy
    d <- dates[min(round, length(dates))]
    sub <- data[data$.date == d, ]

    region_uncertainty <- vapply(regions, function(r) {
      if (has_loc && ".location" %in% names(sub)) {
        rsub <- sub[sub$.location == r, ]
      } else {
        rsub <- sub
      }
      if (nrow(rsub) == 0L) return(1.0)

      freqs <- rsub$.freq
      # Variance of multinomial proportions
      var_est <- sum(freqs * (1 - freqs)) / max(sum(rsub$.count), 1)
      max(var_est, 1e-6)
    }, numeric(1L))

    if (strategy == "thompson") {
      # Thompson: draw from posterior of uncertainty, allocate proportional
      draws <- vapply(seq_along(regions), function(i) {
        stats::rgamma(1, shape = 1 / region_uncertainty[i],
                      rate = cum_alloc[i])
      }, numeric(1L))
      # Invert: higher draw means LESS certain, allocate MORE
      scores <- 1 / (draws + 1e-10)
    } else {
      # UCB: uncertainty + exploration bonus
      scores <- region_uncertainty +
        exploration * sqrt(2 * log(round) / cum_alloc)
    }

    # Allocate proportional to scores
    props <- scores / sum(scores)
    allocs <- as.integer(round(props * capacity))
    # Ensure total = capacity
    diff_total <- capacity - sum(allocs)
    if (diff_total != 0) {
      idx <- which.max(scores)
      allocs[idx] <- allocs[idx] + diff_total
    }
    allocs <- pmax(allocs, 0L)

    cum_alloc <- cum_alloc + allocs

    for (i in seq_along(regions)) {
      freq_est <- if (!is.null(target_lineage)) {
        tl_row <- sub[sub$.lineage == target_lineage, ]
        if (nrow(tl_row) > 0) tl_row$.freq[1] else NA_real_
      } else {
        NA_real_
      }

      alloc_rows <- c(alloc_rows, list(tibble::tibble(
        round       = round,
        region      = regions[i],
        n_allocated = allocs[i],
        uncertainty = region_uncertainty[i],
        frequency   = freq_est
      )))
    }
  }

  alloc_df <- dplyr::bind_rows(alloc_rows)

  summary_df <- alloc_df |>
    dplyr::group_by(.data$region) |>
    dplyr::summarise(
      total_allocated = sum(.data$n_allocated),
      mean_allocated  = mean(.data$n_allocated),
      mean_uncertainty = mean(.data$uncertainty),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$total_allocated))

  structure(
    list(
      allocations = alloc_df,
      summary     = summary_df,
      strategy    = strategy,
      capacity    = capacity
    ),
    class = "adaptive_allocation"
  )
}


#' @export
print.adaptive_allocation <- function(x, ...) {
  cli::cli_h3("Adaptive sequencing allocation")
  cli::cli_text("Strategy: {x$strategy}  |  Capacity: {x$capacity}/round")
  cli::cli_text("{length(unique(x$allocations$round))} rounds, {length(unique(x$allocations$region))} regions")
  cat("\n")
  print(x$summary)
  invisible(x)
}


#' @export
summary.adaptive_allocation <- function(object, ...) {
  object$summary
}


#' Plot adaptive allocation
#' @param x An \code{adaptive_allocation} object.
#' @param ... Unused.
#' @return A ggplot object.
#' @export
plot.adaptive_allocation <- function(x, ...) {
  ggplot2::ggplot(x$allocations,
    ggplot2::aes(x = .data$round, y = .data$n_allocated,
                 fill = .data$region)) +
    ggplot2::geom_col(position = "stack", alpha = 0.8) +
    ggplot2::labs(
      x = "Allocation round", y = "Sequences allocated",
      fill = "Region",
      title = "Adaptive sequencing allocation over time"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}
