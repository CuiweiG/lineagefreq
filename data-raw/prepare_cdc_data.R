# ============================================================
# Prepare CDC SARS-CoV-2 Variant Proportions dataset
#
# Source: CDC COVID Data Tracker
# https://data.cdc.gov/Laboratory-Surveillance/
# SARS-CoV-2-Variant-Proportions/jr58-6ysp
#
# License: Public Domain (U.S. Government Work, 17 USC 105)
#
# Download the full CSV from:
# https://data.cdc.gov/api/views/jr58-6ysp/rows.csv?accessType=DOWNLOAD
#
# Reference:
# Ma KC et al. (2024). Genomic Surveillance for SARS-CoV-2
# Variants. MMWR 73(42). doi:10.15585/mmwr.mm7342a1
# ============================================================

d <- read.csv("SARS-CoV-2_Variant_Proportions.csv",
              stringsAsFactors = FALSE)

# Filter: national level, weighted model, JN.1 emergence period
d <- d[d$usa_or_hhsregion == "USA" &
         d$modeltype == "weighted" &
         d$week_ending >= "2023-10-01" &
         d$week_ending <= "2024-06-30", ]

d$date <- as.Date(substr(d$week_ending, 1, 10))

# Key lineages for the JN.1 emergence narrative
key_variants <- c("JN.1", "HV.1", "EG.5.1", "XBB.1.5", "HK.3",
                  "BA.2.86", "KP.2", "KP.3", "JN.1.11.1")
d$lineage <- ifelse(d$variant %in% key_variants, d$variant, "Other")

# Aggregate shares by date + lineage
agg <- aggregate(share ~ date + lineage, data = d, FUN = sum)

# Convert proportions to approximate counts (8000 seqs/biweek)
total_seqs <- 8000L
result <- do.call(rbind, lapply(split(agg, agg$date), function(x) {
  x$share_norm <- x$share / sum(x$share)
  x$count <- as.integer(round(x$share_norm * total_seqs))
  diff <- total_seqs - sum(x$count)
  if (diff != 0) {
    x$count[which.max(x$count)] <- x$count[which.max(x$count)] + diff
  }
  x
}))

cdc_sarscov2_jn1 <- data.frame(
  date       = result$date,
  lineage    = result$lineage,
  count      = result$count,
  proportion = round(result$share_norm, 6),
  stringsAsFactors = FALSE
)
cdc_sarscov2_jn1 <- cdc_sarscov2_jn1[order(cdc_sarscov2_jn1$date,
                                             cdc_sarscov2_jn1$lineage), ]
rownames(cdc_sarscov2_jn1) <- NULL

usethis::use_data(cdc_sarscov2_jn1, overwrite = TRUE, compress = "xz")
