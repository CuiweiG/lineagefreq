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


#' CDC SARS-CoV-2 variant proportions: JN.1 emergence (US, 2023-2024)
#'
#' A dataset of weekly SARS-CoV-2 variant sequence counts for the
#' United States covering the JN.1 emergence wave (October 2023 to
#' June 2024). Derived from CDC national genomic surveillance
#' proportions, with approximate counts reconstructed using a
#' reference total of 10,000 sequences per week.
#'
#' @format A data frame with 351 rows and 3 columns:
#' \describe{
#'   \item{date}{Week ending date (Date).}
#'   \item{lineage}{Lineage name (character): JN.1, XBB.1.5, EG.5,
#'     HV.1, HK.3, BA.2.86, KP.2, KP.3, Other.}
#'   \item{count}{Approximate sequence count (integer).}
#' }
#'
#' @source Simulated based on published proportions from CDC COVID
#'   Data Tracker and MMWR Vol.73 No.42.
#'   \doi{10.15585/mmwr.mm7342a1}
#'
#' @examples
#' data(cdc_sarscov2_jn1)
#' x <- lfq_data(cdc_sarscov2_jn1,
#'               lineage = lineage, date = date, count = count)
#' x
"cdc_sarscov2_jn1"
