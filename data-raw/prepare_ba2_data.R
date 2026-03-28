# ============================================================
# Prepare CDC BA.1 → BA.2 transition dataset
# Source: CDC COVID Data Tracker (data.cdc.gov)
# License: Public Domain (U.S. Government Work, 17 USC 105)
# Period: US national, Dec 2021 - Jun 2022
# Reference: Lyngse et al. (2022) Lancet Infect Dis
# ============================================================

d <- read.csv("SARS-CoV-2_Variant_Proportions.csv",
              stringsAsFactors = FALSE)
d <- d[d$usa_or_hhsregion == "USA" & d$modeltype == "weighted" &
         d$week_ending >= "2021-12-01" & d$week_ending <= "2022-06-30", ]
d$date <- as.Date(substr(d$week_ending, 1, 10))

key <- c("B.1.1.529", "BA.1", "BA.1.1", "BA.2", "BA.2.12.1", "BA.4", "BA.5")
d$lineage <- ifelse(d$variant %in% key, d$variant, "Other")
d$lineage[d$lineage %in% c("B.1.1.529", "BA.1", "BA.1.1")] <- "BA.1"
d$lineage[d$lineage %in% c("BA.4", "BA.5")] <- "BA.4/5"

agg <- aggregate(share ~ date + lineage, data = d, FUN = sum)
result <- do.call(rbind, lapply(split(agg, agg$date), function(x) {
  x$sn <- x$share / sum(x$share)
  x$count <- as.integer(round(x$sn * 8000))
  diff <- 8000L - sum(x$count)
  if (diff != 0) x$count[which.max(x$count)] <- x$count[which.max(x$count)] + diff
  x
}))

cdc_ba2_transition <- data.frame(
  date = result$date, lineage = result$lineage,
  count = result$count, proportion = round(result$sn, 6),
  stringsAsFactors = FALSE)
cdc_ba2_transition <- cdc_ba2_transition[
  order(cdc_ba2_transition$date, cdc_ba2_transition$lineage), ]
rownames(cdc_ba2_transition) <- NULL
usethis::use_data(cdc_ba2_transition, overwrite = TRUE, compress = "xz")
