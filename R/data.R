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
#' Real surveillance data from the CDC's national genomic surveillance
#' program covering the emergence and dominance of the SARS-CoV-2 JN.1
#' lineage in the United States, October 2023 through June 2024.
#'
#' Derived from CDC's published weighted variant proportion estimates.
#' Approximate biweekly sequence counts were reconstructed from
#' proportions using a reference total of 8,000 sequences per period.
#' The original proportions are retained in the `proportion` column.
#'
#' @format A data frame with 171 rows and 4 columns:
#' \describe{
#'   \item{date}{Biweek ending date (Date).}
#'   \item{lineage}{Lineage name (character): JN.1, XBB.1.5, EG.5.1,
#'     HV.1, HK.3, BA.2.86, KP.2, KP.3, JN.1.11.1, Other.}
#'   \item{count}{Approximate sequence count per biweek (integer).}
#'   \item{proportion}{CDC weighted proportion estimate (numeric).}
#' }
#'
#' @source
#' CDC COVID Data Tracker, SARS-CoV-2 Variant Proportions.
#' Dataset ID: jr58-6ysp.
#' \url{https://data.cdc.gov/Laboratory-Surveillance/SARS-CoV-2-Variant-Proportions/jr58-6ysp}
#'
#' Public domain (U.S. Government Work, 17 USC 105).
#'
#' @references
#' Ma KC, et al. (2024). Genomic Surveillance for SARS-CoV-2 Variants.
#' \emph{MMWR}, 73(42):941--948. \doi{10.15585/mmwr.mm7342a1}
#'
#' @examples
#' data(cdc_sarscov2_jn1)
#' vd <- lfq_data(cdc_sarscov2_jn1,
#'                date = date, lineage = lineage, count = count)
#' fit <- fit_model(vd, engine = "mlr")
#' growth_advantage(fit, type = "relative_Rt", generation_time = 5)
"cdc_sarscov2_jn1"


#' CDC SARS-CoV-2 variant proportions: BA.1 to BA.2 transition (US, 2022)
#'
#' Real surveillance data covering the Omicron BA.1 to BA.2 variant
#' replacement in the United States, December 2021 through June 2022.
#' This is one of the best-documented variant replacement events and
#' serves as an independent validation dataset.
#'
#' @format A data frame with 150 rows and 4 columns:
#' \describe{
#'   \item{date}{Biweek ending date (Date).}
#'   \item{lineage}{Lineage name: BA.1, BA.2, BA.2.12.1, BA.4/5, Other.}
#'   \item{count}{Approximate sequence count per biweek (integer).}
#'   \item{proportion}{CDC weighted proportion estimate (numeric).}
#' }
#'
#' @source CDC COVID Data Tracker (data.cdc.gov, public domain).
#'
#' @references
#' Lyngse FP, et al. (2022). Household transmission of SARS-CoV-2
#' Omicron variant of concern subvariants BA.1 and BA.2 in Denmark.
#' \emph{Nature Communications}, 13:5760.
#' \doi{10.1038/s41467-022-33498-0}
#'
#' @examples
#' data(cdc_ba2_transition)
#' vd <- lfq_data(cdc_ba2_transition,
#'                date = date, lineage = lineage, count = count)
#' fit <- fit_model(vd, engine = "mlr", pivot = "BA.1")
#' # BA.2 Rt ~ 1.34 (consistent with published estimates)
#' growth_advantage(fit, type = "relative_Rt", generation_time = 3.2)
"cdc_ba2_transition"
