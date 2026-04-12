#' Decompose variant fitness into transmissibility and immune escape
#'
#' Partitions the observed growth advantage of each lineage into
#' two components: intrinsic transmissibility (fitness in a fully
#' naive population) and immune escape (additional fitness gained
#' from evading population immunity). This decomposition follows the
#' framework of Figgins and Bedford (2025), where the effective
#' fitness of lineage \eqn{v} is modelled as
#' \eqn{f_v(t) = \beta_v \cdot (1 - \pi_v(t))}
#' with \eqn{\beta_v} representing intrinsic transmissibility and
#' \eqn{\pi_v(t)} the proportion of the population with neutralising
#' immunity against \eqn{v}.
#'
#' @param fit An \code{lfq_fit} object from \code{\link{fit_model}}.
#' @param landscape An \code{immune_landscape} object from
#'   \code{\link{immune_landscape}}.
#' @param generation_time Generation time in days (required for
#'   converting growth rates to reproduction numbers).
#'
#' @return A \code{fitness_decomposition} object (S3 class) with
#'   components:
#'   \describe{
#'     \item{decomposition}{Tibble with columns: \code{lineage},
#'       \code{observed_advantage} (total growth rate), \code{beta}
#'       (intrinsic transmissibility), \code{escape_contribution}
#'       (immune escape component), \code{transmissibility_fraction}
#'       (proportion due to intrinsic fitness),
#'       \code{escape_fraction} (proportion due to immune escape).}
#'     \item{fit}{The input lfq_fit object.}
#'     \item{landscape}{The input immune_landscape.}
#'     \item{generation_time}{Generation time used.}
#'   }
#'
#' @details
#' The decomposition proceeds as follows. For each non-pivot
#' lineage, the observed growth rate \eqn{\delta_v} from the MLR
#' fit is taken as the total fitness difference relative to the
#' reference. The immune escape component is estimated from the
#' immunity differential:
#' \deqn{\text{escape}_v = -\log(1 - \bar{\pi}_v) + \log(1 - \bar{\pi}_{\text{ref}})}
#' where \eqn{\bar{\pi}_v} is the mean population immunity against
#' lineage \eqn{v} over the fitted period. The intrinsic
#' transmissibility component is the residual:
#' \deqn{\beta_v = \delta_v - \text{escape}_v}
#'
#' When cross-immunity data are available in the landscape object,
#' effective immunity against each lineage is computed as the
#' weighted sum of immunity sources, discounted by the cross-immunity
#' matrix.
#'
#' @references
#' Figgins MD, Bedford T (2025). Jointly modeling variant frequencies
#' and case counts to estimate relative variant severity.
#' \emph{medRxiv}. \doi{10.1101/2024.12.02.24318334}
#'
#' @seealso \code{\link{immune_landscape}} for constructing the
#'   immunity input, \code{\link{selective_pressure}} for
#'   population-level early warning signals.
#'
#' @examples
#' \donttest{
#' sim <- simulate_dynamics(n_lineages = 3,
#'   advantages = c("A" = 1.3, "B" = 0.9),
#'   n_timepoints = 15, seed = 1)
#' fit <- fit_model(sim, engine = "mlr")
#'
#' imm_data <- data.frame(
#'   date = rep(unique(sim$.date), each = 3),
#'   lineage = rep(c("A", "B", "ref"), length(unique(sim$.date))),
#'   immunity = c(rep(c(0.2, 0.5, 0.4),
#'     length(unique(sim$.date))))
#' )
#' il <- immune_landscape(imm_data, date, lineage, immunity)
#' fd <- fitness_decomposition(fit, il, generation_time = 5)
#' fd
#' }
#'
#' @export
fitness_decomposition <- function(fit, landscape, generation_time) {

  if (!inherits(fit, "lfq_fit")) {
    cli::cli_abort("{.arg fit} must be an {.cls lfq_fit} object.")
  }
  if (!inherits(landscape, "immune_landscape")) {
    cli::cli_abort("{.arg landscape} must be an {.cls immune_landscape} object.")
  }
  if (missing(generation_time) || !is.numeric(generation_time) ||
      generation_time <= 0) {
    cli::cli_abort("{.arg generation_time} must be a positive number.")
  }

  pivot    <- fit$pivot
  lineages <- fit$lineages
  non_pivot <- setdiff(lineages, pivot)
  ts       <- fit$time_scale
  tau      <- generation_time / ts

  # Get growth rates
  deltas <- fit$growth_rates

  # Compute mean immunity per lineage over the fitted period
  imm_est <- landscape$estimates

  # Match lineages between fit and landscape
  common <- intersect(lineages, landscape$lineages)
  if (length(common) < 2L) {
    cli::cli_abort(
      "Need at least 2 lineages in common between fit and landscape; found {length(common)}."
    )
  }

  mean_immunity <- vapply(lineages, function(lin) {
    sub <- imm_est[imm_est$lineage == lin, ]
    if (nrow(sub) == 0L) return(0)
    mean(sub$immunity, na.rm = TRUE)
  }, numeric(1L))

  # Apply cross-immunity if available
  if (!is.null(landscape$cross_immunity)) {
    xim <- landscape$cross_immunity
    # Effective immunity: sum over source lineages weighted by cross-immunity
    eff_immunity <- vapply(lineages, function(target) {
      if (!target %in% rownames(xim)) return(mean_immunity[target])
      sources <- intersect(names(mean_immunity), colnames(xim))
      weights <- xim[sources, target]
      sum(mean_immunity[sources] * weights)
    }, numeric(1L))
    mean_immunity <- pmin(eff_immunity, 1)
  }

  # Decomposition
  # Effective susceptibility: S_v = 1 - pi_v
  # Escape component: differential susceptibility in log space
  # escape_v = -log(S_v) + log(S_ref) = log(S_ref / S_v)
  s_ref <- 1 - mean_immunity[pivot]
  s_ref <- max(s_ref, 0.01)  # floor to avoid log(0)

  rows <- list()
  for (lin in non_pivot) {
    delta <- deltas[lin]
    s_v <- max(1 - mean_immunity[lin], 0.01)

    escape_component <- log(s_ref / s_v) * tau
    intrinsic_component <- delta - escape_component

    total_abs <- abs(intrinsic_component) + abs(escape_component)
    if (total_abs < 1e-10) {
      trans_frac <- 0.5
      esc_frac   <- 0.5
    } else {
      trans_frac <- abs(intrinsic_component) / total_abs
      esc_frac   <- abs(escape_component) / total_abs
    }

    rows <- c(rows, list(tibble::tibble(
      lineage                  = lin,
      observed_advantage       = delta,
      beta                     = intrinsic_component,
      escape_contribution      = escape_component,
      transmissibility_fraction = trans_frac,
      escape_fraction          = esc_frac
    )))
  }

  # Add pivot row
  rows <- c(list(tibble::tibble(
    lineage = pivot,
    observed_advantage = 0,
    beta = 0,
    escape_contribution = 0,
    transmissibility_fraction = NA_real_,
    escape_fraction = NA_real_
  )), rows)

  decomp <- dplyr::bind_rows(rows)


  structure(
    list(
      decomposition   = decomp,
      fit             = fit,
      landscape       = landscape,
      generation_time = generation_time
    ),
    class = "fitness_decomposition"
  )
}


#' @export
print.fitness_decomposition <- function(x, ...) {
  cli::cli_h3("Fitness decomposition")
  cli::cli_text(
    "Pivot: {x$fit$pivot}  |  Generation time: {x$generation_time} days"
  )
  cat("\n")
  d <- x$decomposition[!is.na(x$decomposition$transmissibility_fraction), ]
  for (i in seq_len(nrow(d))) {
    row <- d[i, ]
    cli::cli_text(
      "{row$lineage}: {round(row$transmissibility_fraction * 100)}% intrinsic / {round(row$escape_fraction * 100)}% escape (delta = {round(row$observed_advantage, 4)})"
    )
  }
  invisible(x)
}


#' Tidy a fitness decomposition
#'
#' @param x A \code{fitness_decomposition} object.
#' @param ... Unused.
#' @return A tibble with the decomposition results.
#' @export
tidy.fitness_decomposition <- function(x, ...) {
  x$decomposition
}


#' Plot fitness decomposition
#'
#' Stacked bar plot showing the partition of growth advantage into
#' intrinsic transmissibility and immune escape components.
#'
#' @param x A \code{fitness_decomposition} object.
#' @param ... Unused.
#' @return A ggplot object.
#' @export
plot.fitness_decomposition <- function(x, ...) {

  d <- x$decomposition[x$decomposition$lineage != x$fit$pivot, ]

  plot_data <- tidyr::pivot_longer(
    d[, c("lineage", "beta", "escape_contribution")],
    cols = c("beta", "escape_contribution"),
    names_to = "component",
    values_to = "value"
  )

  plot_data$component <- ifelse(
    plot_data$component == "beta",
    "Intrinsic transmissibility",
    "Immune escape"
  )

  ggplot2::ggplot(plot_data,
    ggplot2::aes(x = .data$lineage, y = .data$value,
                 fill = .data$component)) +
    ggplot2::geom_col(position = "stack", alpha = 0.85) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                        colour = "grey30") +
    ggplot2::labs(
      x = "Lineage", y = "Growth rate component",
      fill = "Component",
      title = "Fitness decomposition"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}
