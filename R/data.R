#' Simulated SARS-CoV-2 variant frequency data (US, 2022)
#'
#' A simulated dataset of weekly SARS-CoV-2 variant sequence counts
#' for the United States in 2022. Includes the BA.1 to BA.2 to
#' BA.4/5 to BQ.1 transition dynamics. This is simulated data (not
#' real GISAID data) to avoid license restrictions while preserving
#' realistic statistical properties.
#'
#' @format A data frame with 200 rows and 4 columns:
#' \describe{
#'   \item{date}{Collection date (Date, weekly).}
#'   \item{variant}{Variant name (character): BA.1, BA.2, BA.4/5,
#'     BQ.1, Other.}
#'   \item{count}{Number of sequences assigned to this variant
#'     in this week (integer).}
#'   \item{total}{Total sequences for this week (integer).}
#' }
#'
#' @source Simulated based on parameters from published CDC MMWR
#'   genomic surveillance reports and 'Nextstrain' public data.
#'
#' @examples
#' data(sarscov2_us_2022)
#' x <- lfq_data(sarscov2_us_2022, lineage = variant,
#'               date = date, count = count, total = total)
#' x
"sarscov2_us_2022"


#' Simulated influenza A/H3N2 clade frequency data
#'
#' A simulated dataset of weekly influenza A/H3N2 clade sequence
#' counts over a single Northern Hemisphere season (24 weeks).
#'
#' @format A data frame with 96 rows and 4 columns:
#' \describe{
#'   \item{date}{Collection date (Date, weekly).}
#'   \item{clade}{Clade name (character).}
#'   \item{count}{Number of sequences (integer).}
#'   \item{total}{Total sequences for this week (integer).}
#' }
#'
#' @source Simulated based on 'Nextstrain' influenza clade dynamics.
#'
#' @examples
#' data(influenza_h3n2)
#' x <- lfq_data(influenza_h3n2, lineage = clade,
#'               date = date, count = count, total = total)
#' x
"influenza_h3n2"
