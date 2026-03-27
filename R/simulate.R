#' Simulate lineage frequency dynamics
#'
#' Generates synthetic lineage frequency data under a multinomial
#' sampling model with configurable growth advantages. Useful for
#' model validation, power analysis, and teaching.
#'
#' @param n_lineages Number of lineages including reference. Default 4.
#' @param n_timepoints Number of time points. Default 20.
#' @param total_per_tp Sequences per time point. A single integer
#'   (constant across time) or a vector of length `n_timepoints`
#'   (variable sampling effort). Default 500.
#' @param advantages Named numeric vector of per-week multiplicative
#'   growth advantages for non-reference lineages. Length must equal
#'   `n_lineages - 1`. Values > 1 mean growing faster than reference,
#'   < 1 means declining. If `NULL`, random values in (0.8, 1.5).
#' @param reference_name Name of the reference lineage. Default `"ref"`.
#' @param start_date Start date for the time series. Default today.
#' @param interval Days between consecutive time points. Default 7
#'   (weekly data).
#' @param overdispersion If `NULL` (default), standard multinomial
#'   sampling. If a positive number, Dirichlet-Multinomial sampling
#'   with concentration = 1 / overdispersion. Larger values = more
#'   overdispersion.
#' @param seed Random seed for reproducibility.
#'
#' @return An [lfq_data] object with an additional `true_freq` column
#'   containing the true (pre-sampling) frequencies.
#'
#' @examples
#' # JN.1 grows 1.3x per week, KP.3 declines at 0.9x
#' sim <- simulate_dynamics(
#'   n_lineages = 3,
#'   advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
#'   n_timepoints = 15,
#'   seed = 42
#' )
#' sim
#'
#' @export
simulate_dynamics <- function(n_lineages    = 4L,
                              n_timepoints  = 20L,
                              total_per_tp  = 500L,
                              advantages    = NULL,
                              reference_name = "ref",
                              start_date    = Sys.Date(),
                              interval      = 7L,
                              overdispersion = NULL,
                              seed          = NULL) {

  assert_pos_int(n_lineages,   "n_lineages")
  assert_pos_int(n_timepoints, "n_timepoints")

  if (!is.null(seed)) set.seed(seed)

  # --- Totals per time point ---
  if (length(total_per_tp) == 1L) {
    totals <- rep(as.integer(total_per_tp), n_timepoints)
  } else {
    if (length(total_per_tp) != n_timepoints) {
      cli::cli_abort(
        "{.arg total_per_tp} must be length 1 or {n_timepoints}."
      )
    }
    totals <- as.integer(total_per_tp)
  }

  # --- Growth advantages ---
  if (is.null(advantages)) {
    advantages <- stats::runif(n_lineages - 1L, min = 0.8, max = 1.5)
  }
  if (length(advantages) != n_lineages - 1L) {
    cli::cli_abort(
      "{.arg advantages} must have length {n_lineages - 1} (one per non-reference lineage)."
    )
  }
  if (is.null(names(advantages))) {
    names(advantages) <- paste0("V", seq_along(advantages) + 1L)
  }

  # --- Lineage names and per-step log growth rates ---
  lin_names <- c(reference_name, names(advantages))
  deltas    <- c(0, log(advantages) / 7 * interval)

  # --- Initial log-proportions (equal) ---
  init_lp <- rep(0, n_lineages)

  # --- Generate dates ---
  dates <- start_date + (seq_len(n_timepoints) - 1L) * interval

  # --- Simulate ---
  freq_mat  <- matrix(NA_real_,    nrow = n_timepoints, ncol = n_lineages)
  count_mat <- matrix(NA_integer_, nrow = n_timepoints, ncol = n_lineages)

  for (t in seq_len(n_timepoints)) {
    lp    <- init_lp + deltas * (t - 1L)
    probs <- exp(lp - log_sum_exp(lp))
    freq_mat[t, ] <- probs

    # Ensure strictly positive + normalised
    probs <- pmax(probs, .Machine$double.eps)
    probs <- probs / sum(probs)

    if (is.null(overdispersion)) {
      count_mat[t, ] <- stats::rmultinom(1L, size = totals[t], prob = probs)
    } else {
      # Dirichlet-Multinomial via Gamma mixture
      alpha <- probs / overdispersion
      g     <- stats::rgamma(n_lineages, shape = alpha, rate = 1)
      g     <- g / sum(g)
      count_mat[t, ] <- stats::rmultinom(1L, size = totals[t], prob = g)
    }
  }

  # --- Assemble data frame ---
  df <- data.frame(
    date      = rep(dates, each = n_lineages),
    lineage   = rep(lin_names, n_timepoints),
    count     = as.vector(t(count_mat)),
    true_freq = as.vector(t(freq_mat)),
    stringsAsFactors = FALSE
  )

  lfq_data(df, lineage = lineage, date = date, count = count)
}
