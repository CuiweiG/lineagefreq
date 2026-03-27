#' @title lineagefreq: Lineage Frequency Dynamics from Genomic
#' Surveillance Counts
#'
#' @description
#' Models pathogen lineage frequency dynamics from genomic
#' surveillance count data. Provides a unified `fit_model()`
#' interface with multiple engines (MLR, hierarchical MLR,
#' Piantham), rolling-origin backtesting via `backtest()`,
#' standardized scoring via `score_forecasts()`, emergence
#' detection, and sequencing power analysis.
#'
#' @section Quick start:
#' ```
#' x <- lfq_data(my_counts, lineage = lineage, date = date, count = n)
#' fit <- fit_model(x, engine = "mlr")
#' growth_advantage(fit, generation_time = 5)
#' ```
#'
#' @importFrom rlang .data %||%
#' @importFrom stats optim setNames coef vcov logLik nobs AIC
#'   confint predict residuals rmultinom dmultinom qnorm median
#'   pchisq quantile p.adjust glm binomial rgamma dnorm
#' @importFrom ggplot2 autoplot
#' @keywords internal
"_PACKAGE"

#' @return See the specific method documentation.
#' @export
ggplot2::autoplot

#' Package version and system information
#'
#' Reports lineagefreq version and availability of optional
#' backends. Useful for reproducibility and bug reports.
#'
#' @return A list with components: `version`, `r_version`,
#'   `stan_available`, `engines`.
#'
#' @examples
#' lfq_version()
#'
#' @export
lfq_version <- function() {
  eng <- lfq_engines()
  list(
    version        = as.character(utils::packageVersion("lineagefreq")),
    r_version      = paste0(R.version$major, ".", R.version$minor),
    stan_available = lfq_stan_available(),
    engines        = eng$engine[eng$available]
  )
}

# False-positive NSE bindings from data.frame column references
utils::globalVariables(c("lineage", "count", "date", "variant",
                         "clade", "total", "location", "t_num"))
