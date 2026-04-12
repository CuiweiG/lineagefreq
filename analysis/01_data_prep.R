###############################################################################
# 01_data_prep.R — All data preprocessing
# lineagefreq validation analysis
###############################################################################

cat("[01_data_prep] Starting...\n")
source("analysis/00_setup.R")

if (packageVersion("lineagefreq") < "0.5.0") {
  stop("Please install local dev version: devtools::install('C:/Users/cg223/Desktop/lineagefreq')")
}

###############################################################################
# Section A: ECDC multi-country data
###############################################################################

cat("  [A] ECDC multi-country data...\n")

ecdc_raw <- read_csv("analysis/data/ecdc_variants.csv", show_col_types = FALSE)
cat(sprintf("    Raw ECDC: %d rows, %d columns\n", nrow(ecdc_raw), ncol(ecdc_raw)))

# Five countries with best genomic surveillance coverage in Europe
target_countries <- c("Denmark", "Germany", "France", "Netherlands", "Spain")

ecdc_filtered <- ecdc_raw |>
  filter(country %in% target_countries)

cat(sprintf("    After country filter: %d rows\n", nrow(ecdc_filtered)))

# Convert year_week to Date (Monday of ISO week)
# ECDC format is "2022-01" (no W prefix), but ISOweek2date requires "2022-W01-1"
cat(sprintf("    Sample year_week values: %s\n",
            paste(head(ecdc_filtered$year_week, 3), collapse = ", ")))
ecdc_filtered <- ecdc_filtered |>
  mutate(
    date = ISOweek::ISOweek2date(paste0(sub("-", "-W", year_week), "-1"))
  )

# Define analysis periods
# BA.2 transition: Dec 2021 – Jun 2022 (BA.1 → BA.2 replacement)
# JN.1 emergence:  Oct 2023 – Mar 2024 (JN.1 rapid expansion)
ba2_start <- as.Date("2021-12-01")
ba2_end   <- as.Date("2022-06-30")
jn1_start <- as.Date("2023-10-01")
jn1_end   <- as.Date("2024-03-31")

ecdc_ba2 <- ecdc_filtered |>
  filter(date >= ba2_start, date <= ba2_end)

ecdc_jn1 <- ecdc_filtered |>
  filter(date >= jn1_start, date <= jn1_end)

cat(sprintf("    BA.2 period: %d rows (%s to %s)\n",
            nrow(ecdc_ba2), ba2_start, ba2_end))
cat(sprintf("    JN.1 period: %d rows (%s to %s)\n",
            nrow(ecdc_jn1), jn1_start, jn1_end))

# Rename columns to lineagefreq conventions
rename_ecdc <- function(df) {
  df |>
    rename(
      lineage = variant,
      count   = number_detections_variant
    ) |>
    # Ensure required columns exist; adapt to actual ECDC column names
    # If column names differ, try alternative names
    mutate(
      count = as.numeric(count)
    ) |>
    filter(!is.na(count), count >= 0)
}

ecdc_ba2 <- rename_ecdc(ecdc_ba2)
ecdc_jn1 <- rename_ecdc(ecdc_jn1)

# Create lfq_data objects per country-period
create_lfq_objects <- function(df, period_label) {
  countries <- unique(df$country)
  result <- list()

  for (cntry in countries) {
    cntry_data <- df |> filter(country == cntry)

    tryCatch({
      # Aggregate by date and lineage; keep long format for lfq_data()
      agg <- cntry_data |>
        group_by(date, lineage) |>
        summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

      obj <- lfq_data(agg, lineage = lineage, date = date, count = count)

      key <- paste0(cntry, "_", period_label)
      result[[key]] <- obj
      cat(sprintf("    Created lfq_data: %s (%d dates, %d lineages)\n",
                  key, n_distinct(agg$date), n_distinct(agg$lineage)))
    },
    error = function(e) {
      cat(sprintf("    WARNING: Failed to create lfq_data for %s_%s: %s\n",
                  cntry, period_label, conditionMessage(e)))
    })
  }
  result
}

ecdc_objects_ba2 <- create_lfq_objects(ecdc_ba2, "BA2")
ecdc_objects_jn1 <- create_lfq_objects(ecdc_jn1, "JN1")
ecdc_prepared    <- c(ecdc_objects_ba2, ecdc_objects_jn1)

saveRDS(ecdc_prepared, "analysis/results/ecdc_prepared.rds")
cat(sprintf("    Saved %d lfq_data objects to ecdc_prepared.rds\n",
            length(ecdc_prepared)))

###############################################################################
# Section B: CDC HHS regional data
###############################################################################

cat("  [B] CDC HHS regional data...\n")

cdc_raw <- read_csv("analysis/data/cdc_variant_proportions_full.csv",
                     show_col_types = FALSE)
cat(sprintf("    Raw CDC: %d rows, %d columns\n", nrow(cdc_raw), ncol(cdc_raw)))

# Filter to weighted model and biweekly intervals
cdc_filtered <- cdc_raw |>
  filter(
    modeltype == "weighted",
    grepl("biweek", time_interval, ignore.case = TRUE)
  )

# Extract HHS region number; keep regions 1-10
cdc_filtered <- cdc_filtered |>
  mutate(
    region_num = as.integer(str_extract(usa_or_hhsregion, "\\d+"))
  ) |>
  filter(!is.na(region_num), region_num >= 1, region_num <= 10)

cat(sprintf("    After filtering: %d rows across %d regions\n",
            nrow(cdc_filtered), n_distinct(cdc_filtered$region_num)))

# BA.2 period: Jan 2022 – Jun 2022
cdc_ba2 <- cdc_filtered |>
  filter(
    week_ending >= as.Date("2022-01-01"),
    week_ending <= as.Date("2022-06-30")
  )

# Reconstruct counts from proportions
# ASSUMPTION: ~800 sequences per HHS region per biweekly period.
# Justification: CDC Genomic Surveillance dashboard reports ~50,000-100,000
# sequences nationally per month during Omicron surge; 10 regions × 2 biweeks
# ≈ ~2500-5000 per region-biweek. We use 800 as a conservative lower bound
# to avoid overconfidence in count-based models. This is a known limitation.
ASSUMED_SEQS_PER_BIWEEK <- 800L

cdc_regional <- list()

for (region in 1:10) {
  region_data <- cdc_ba2 |>
    filter(region_num == region)

  if (nrow(region_data) == 0) {
    cat(sprintf("    WARNING: No data for Region %d\n", region))
    next
  }

  tryCatch({
    # Reconstruct counts from proportions in long format
    region_long <- region_data |>
      select(date = week_ending, lineage = variant, share) |>
      mutate(
        date  = as.Date(date),
        share = replace_na(share, 0),
        count = as.integer(round(share * ASSUMED_SEQS_PER_BIWEEK))
      ) |>
      select(date, lineage, count)

    obj <- lfq_data(region_long, lineage = lineage, date = date, count = count)
    key <- paste0("Region_", region)
    cdc_regional[[key]] <- obj
    cat(sprintf("    Region %d: %d dates, %d lineages\n",
                region, n_distinct(region_long$date),
                n_distinct(region_long$lineage)))
  },
  error = function(e) {
    cat(sprintf("    WARNING: Failed for Region %d: %s\n",
                region, conditionMessage(e)))
  })
}

saveRDS(cdc_regional, "analysis/results/cdc_regional_prepared.rds")
cat(sprintf("    Saved %d regional lfq_data objects\n", length(cdc_regional)))

###############################################################################
# Section C: Immune landscape construction
###############################################################################

cat("  [C] Immune landscape construction...\n")

# ── C.1: Natural infection proxy from Anti-N seroprevalence ──────────────────

sero_raw <- read_csv("analysis/data/cdc_seroprevalence.csv",
                      show_col_types = FALSE)
cat(sprintf("    Seroprevalence: %d rows\n", nrow(sero_raw)))

# Extract national Anti-N All Ages seroprevalence by round/date
# Actual CSV columns include:
#   Site, Date Range of Specimen Collection, Round,
#   "Rate (%) [Anti-N, All Ages Cumulative Prevalence, Rounds 1-30 only]"
# There is no age_group column — the age group is encoded in column names.

# Identify the Anti-N All Ages rate column
anti_n_all_ages_col <- grep(
  "Rate.*Anti-N.*All Ages", names(sero_raw), value = TRUE
)
if (length(anti_n_all_ages_col) == 0) {
  # Fallback: try any Anti-N rate column

  anti_n_all_ages_col <- grep("Rate.*Anti-N", names(sero_raw), value = TRUE)
}
cat(sprintf("    Using sero column: %s\n", anti_n_all_ages_col[1]))

# ── Diagnostics: inspect CSV structure ────────────────────────────────────────
cat(sprintf("    Unique Site values: %s\n",
            paste(head(unique(sero_raw$Site), 20), collapse = ", ")))
cat(sprintf("    Columns containing 'Anti' or 'Rate': %s\n",
            paste(grep("Anti|Rate|rate", names(sero_raw), value = TRUE),
                  collapse = ", ")))

# Parse the date range: extract the end date from "Mon DD - Mon DD, YYYY"
# e.g. "Aug 6 - Aug 11, 2020" → "Aug 11, 2020"
# Site column uses state abbreviations + "US" for national
sero_national <- sero_raw |>
  filter(Site == "US") |>
  mutate(
    date_str = sub("^.*-\\s*", "", `Date Range of Specimen Collection`),
    date = as.Date(date_str, format = "%b %d, %Y"),
    # Rate column may contain non-numeric characters; use parse_number()
    seroprevalence_pct = readr::parse_number(
      as.character(.data[[anti_n_all_ages_col[1]]])
    ),
    natural_infection_rate = seroprevalence_pct / 100
  ) |>
  filter(!is.na(date), !is.na(natural_infection_rate)) |>
  select(date, seroprevalence_pct, natural_infection_rate) |>
  arrange(date)

cat(sprintf("    Filtered rows: %d\n", nrow(sero_national)))
cat(sprintf("    natural_infection_rate values: %s\n",
            paste(head(sero_national$natural_infection_rate, 10), collapse = ", ")))
cat(sprintf("    Date values: %s\n",
            paste(head(sero_national$date, 10), collapse = ", ")))

# ── Fallback: if CSV parsing yielded no usable rows, use published estimates ──
# Sources: Clarke et al. 2022 JAMA; Jones et al. 2023 MMWR
if (nrow(sero_national) == 0 ||
    all(is.na(sero_national$natural_infection_rate))) {
  cat("    WARNING: CSV parsing returned 0 usable rows.\n")
  cat("    Using published CDC seroprevalence estimates (Anti-N, All Ages).\n")

  sero_national <- tibble(
    date = as.Date(c(
      "2020-07-15",   # Round 1
      "2020-10-15",   # Round 4
      "2021-06-15",   # Round 10
      "2021-12-15",   # Round 15
      "2022-02-15",   # Round 18
      "2022-05-15",   # Round 21
      "2022-08-15"    # Round 24
    )),
    seroprevalence_pct = c(3.5, 7.0, 20.0, 34.0, 58.0, 73.0, 79.0),
    natural_infection_rate = c(3.5, 7.0, 20.0, 34.0, 58.0, 73.0, 79.0) / 100
  )
}

cat(sprintf("    National sero rounds: %d rows, date range: %s to %s\n",
            nrow(sero_national),
            min(sero_national$date), max(sero_national$date)))
cat(sprintf("    Seroprevalence range: %.1f%% to %.1f%%\n",
            min(sero_national$natural_infection_rate, na.rm = TRUE) * 100,
            max(sero_national$natural_infection_rate, na.rm = TRUE) * 100))

# ── C.2: Vaccination data from OWID ─────────────────────────────────────────

vacc_raw <- read_csv("analysis/data/owid_vaccinations.csv",
                      show_col_types = FALSE)

vacc_us <- vacc_raw |>
  filter(location == "United States") |>
  select(date, people_fully_vaccinated_per_hundred,
         total_boosters_per_hundred) |>
  mutate(date = as.Date(date)) |>
  filter(!is.na(date)) |>
  arrange(date) |>
  # Forward-fill missing values

  fill(people_fully_vaccinated_per_hundred, total_boosters_per_hundred,
       .direction = "down") |>
  mutate(
    people_fully_vaccinated_per_hundred = replace_na(people_fully_vaccinated_per_hundred, 0),
    total_boosters_per_hundred          = replace_na(total_boosters_per_hundred, 0)
  )

cat(sprintf("    US vaccination data: %d days\n", nrow(vacc_us)))

# ── C.3: Construct time-varying immunity for BA.2 period ─────────────────────
# Reference period: Jan 2022 – Jun 2022
# Vaccine efficacy waning model:
#   - Half-life of vaccine protection: 6 months (Goldberg et al. 2022, NEJM)
#   - Peak US vaccination: Mar-May 2021 (2-dose), boosters from Sep 2021
#   - By Jan 2022 (~10 months post peak 2-dose): efficacy ~35%
#   - By Jan 2022 (~3 months post booster start): booster efficacy ~75%
# Variant-specific immune escape:
#   - BA.2: ~20% escape (Tuekprakhon et al. 2022, Nature)
#   - BA.4/5: ~40% escape (Tuekprakhon et al. 2022, Nature)

immunity_dates <- seq(as.Date("2022-01-01"), as.Date("2022-06-30"), by = "week")

# Interpolate seroprevalence to weekly
sero_interp <- approx(
  x = as.numeric(sero_national$date),
  y = sero_national$natural_infection_rate,
  xout = as.numeric(immunity_dates),
  rule = 2  # constant extrapolation outside range
)$y

# Get vaccination rates at each date
vacc_at_dates <- vacc_us |>
  filter(date %in% immunity_dates | date >= min(immunity_dates))

vacc_interp_full <- approx(
  x = as.numeric(vacc_us$date),
  y = vacc_us$people_fully_vaccinated_per_hundred / 100,
  xout = as.numeric(immunity_dates),
  rule = 2
)$y

vacc_interp_boost <- approx(
  x = as.numeric(vacc_us$date),
  y = vacc_us$total_boosters_per_hundred / 100,
  xout = as.numeric(immunity_dates),
  rule = 2
)$y

# Waning model: exponential decay with 6-month half-life
# Reference vaccination peak: 2021-04-15 (2-dose), 2021-11-01 (booster)
waning_halflife_days <- 180
peak_2dose  <- as.Date("2021-04-15")
peak_boost  <- as.Date("2021-11-01")

vaccine_eff <- numeric(length(immunity_dates))
for (i in seq_along(immunity_dates)) {
  days_since_2dose  <- as.numeric(immunity_dates[i] - peak_2dose)
  days_since_boost  <- as.numeric(immunity_dates[i] - peak_boost)

  eff_2dose  <- 0.90 * 2^(-days_since_2dose / waning_halflife_days)
  eff_boost  <- 0.95 * 2^(-max(0, days_since_boost) / waning_halflife_days)

  # Weighted by fraction boosted vs 2-dose only
  frac_boosted <- min(1, vacc_interp_boost[i] / max(vacc_interp_full[i], 0.01))
  frac_2dose   <- 1 - frac_boosted

  vaccine_eff[i] <- vacc_interp_full[i] * (frac_2dose * eff_2dose +
                                             frac_boosted * eff_boost)
}

# Total effective immunity: 1 - (1 - vaccine_eff) * (1 - infection_eff)
total_immunity <- 1 - (1 - vaccine_eff) * (1 - sero_interp)

# Variant-specific escape factors
# BA.2: ~20% escape → effective immunity reduced by 20%
# BA.4/5: ~40% escape → effective immunity reduced by 40%
# Reference: Tuekprakhon et al. 2022, Nature
escape_ba2  <- 0.20
escape_ba45 <- 0.40

immunity_vs_ba2  <- total_immunity * (1 - escape_ba2)
immunity_vs_ba45 <- total_immunity * (1 - escape_ba45)

immune_df <- tibble(
  date              = immunity_dates,
  natural_infection  = sero_interp,
  vaccine_effective  = vaccine_eff,
  total_immunity     = total_immunity,
  immunity_vs_ba2    = immunity_vs_ba2,
  immunity_vs_ba45   = immunity_vs_ba45
)

cat(sprintf("    Immunity range: %.1f%% to %.1f%% (total, Jan-Jun 2022)\n",
            min(total_immunity) * 100, max(total_immunity) * 100))

# Create immune_landscape object
# API: immune_landscape(data, date, lineage, immunity, type, cross_immunity)
# Expects a long-format data frame with date, lineage, and immunity columns.
# We create one row per date × variant with variant-specific effective immunity.
immune_long <- bind_rows(
  tibble(date = immunity_dates, lineage = "BA.2",  immunity = immunity_vs_ba2),
  tibble(date = immunity_dates, lineage = "BA.4/5", immunity = immunity_vs_ba45)
)

tryCatch({
  immune_landscape_obj <- immune_landscape(
    data     = immune_long,
    date     = date,
    lineage  = lineage,
    immunity = immunity
  )
  cat("    Created immune_landscape object\n")
},
error = function(e) {
  cat(sprintf("    WARNING: immune_landscape() failed: %s\n",
              conditionMessage(e)))
  cat("    Saving raw data frame instead; 05_fitness.R will construct object\n")
  immune_landscape_obj <<- immune_df
})

saveRDS(immune_landscape_obj, "analysis/results/immune_landscape.rds")
cat("    Saved immune_landscape.rds\n")

###############################################################################
# Section D: WHO FluNet data
###############################################################################

cat("  [D] WHO FluNet data...\n")

# FluNet data format can vary. Attempt to read as CSV; if HTML or incompatible,
# fall back to package built-in simulated data.

flunet_prepared <- list()

tryCatch({
  flu_raw <- read_csv("analysis/data/who_flunet.csv", show_col_types = FALSE)
  cat(sprintf("    Raw FluNet: %d rows, %d columns\n", nrow(flu_raw), ncol(flu_raw)))

  # FluNet columns: COUNTRY_AREA_TERRITORY, ISO_YEAR, ISO_WEEK,
  # AH1N12009, AH3, INF_B, BVIC_2DEL, BYAM, etc.
  flu_cols <- names(flu_raw)
  cat(sprintf("    Columns: %s\n", paste(head(flu_cols, 15), collapse = ", ")))

  # Print available countries for diagnostics
  cat(sprintf("    Sample countries: %s\n",
              paste(head(unique(flu_raw$COUNTRY_AREA_TERRITORY), 10), collapse = ", ")))

  # Target countries with good surveillance
  flu_countries <- c("United States of America", "Australia",
                     "United Kingdom of Great Britain and Northern Ireland")

  # Filter to available countries and recent seasons (2022-2024)
  flu_filtered <- flu_raw |>
    filter(
      COUNTRY_AREA_TERRITORY %in% flu_countries,
      ISO_YEAR >= 2022
    )

  if (nrow(flu_filtered) == 0) {
    # Try partial match if exact names differ
    flu_filtered <- flu_raw |>
      filter(grepl("United States|Australia|United Kingdom",
                   COUNTRY_AREA_TERRITORY, ignore.case = TRUE),
             ISO_YEAR >= 2022)
  }

  if (nrow(flu_filtered) > 0) {
    cat(sprintf("    Filtered FluNet: %d rows\n", nrow(flu_filtered)))

    # Extract subtype counts; actual FluNet columns:
    # AH1N12009, AH3, INF_B, BVIC_2DEL, BVIC_3DEL, BYAM
    # Use specific subtypes for cleaner analysis
    subtype_cols <- intersect(
      c("AH1N12009", "AH3", "BVIC_2DEL", "BYAM", "INF_B"),
      flu_cols
    )
    cat(sprintf("    Subtype columns found: %s\n",
                paste(subtype_cols, collapse = ", ")))

    if (length(subtype_cols) >= 2) {
      for (cntry in unique(flu_filtered$COUNTRY_AREA_TERRITORY)) {
        cntry_data <- flu_filtered |>
          filter(COUNTRY_AREA_TERRITORY == cntry)

        # Construct date from ISO_YEAR + ISO_WEEK
        cntry_data <- cntry_data |>
          mutate(
            date = ISOweek::ISOweek2date(
              sprintf("%04d-W%02d-1", as.integer(ISO_YEAR), as.integer(ISO_WEEK))
            )
          )

        # Pivot subtype counts to long format for lfq_data()
        flu_long <- cntry_data |>
          select(date, all_of(subtype_cols)) |>
          pivot_longer(cols = all_of(subtype_cols),
                       names_to = "lineage", values_to = "count") |>
          mutate(count = replace_na(as.integer(as.numeric(count)), 0L)) |>
          filter(count >= 0)

        tryCatch({
          obj <- lfq_data(flu_long, lineage = lineage, date = date, count = count)
          cntry_clean <- gsub(" ", "_", sub(" of America", "", cntry))
          flunet_prepared[[paste0("flu_", cntry_clean)]] <- obj
          cat(sprintf("    Created flu lfq_data: %s (%d weeks, %d subtypes)\n",
                      cntry_clean, n_distinct(flu_long$date),
                      n_distinct(flu_long$lineage)))
        },
        error = function(e) {
          cat(sprintf("    WARNING: lfq_data failed for %s: %s\n",
                      cntry, conditionMessage(e)))
        })
      }
    } else {
      cat("    WARNING: Could not identify subtype columns in FluNet data\n")
      cat("    Available columns: ", paste(flu_cols, collapse = ", "), "\n")
    }
  } else {
    cat("    WARNING: No matching countries found in FluNet data\n")
  }
},
error = function(e) {
  cat(sprintf("    FluNet read failed: %s\n", conditionMessage(e)))
})

# Fallback: if no real flu data, note for downstream scripts
if (length(flunet_prepared) == 0) {
  cat("    NOTE: Real FluNet data unavailable or incompatible.\n")
  cat("    07_influenza.R will use package built-in simulated data.\n")
  attr(flunet_prepared, "source") <- "simulated"
} else {
  attr(flunet_prepared, "source") <- "who_flunet"
}

saveRDS(flunet_prepared, "analysis/results/flunet_prepared.rds")
cat(sprintf("    Saved flunet_prepared.rds (%d datasets)\n",
            length(flunet_prepared)))

###############################################################################
# Summary
###############################################################################

cat("\n  ══════════════════════════════════════════\n")
cat("  Data preparation summary:\n")
cat(sprintf("    ECDC:     %d country-period datasets\n", length(ecdc_prepared)))
cat(sprintf("    CDC HHS:  %d regional datasets\n", length(cdc_regional)))
cat(sprintf("    Immunity: %d weekly time points (%.1f%% → %.1f%%)\n",
            nrow(immune_df),
            min(immune_df$total_immunity) * 100,
            max(immune_df$total_immunity) * 100))
cat(sprintf("    FluNet:   %d datasets (source: %s)\n",
            length(flunet_prepared), attr(flunet_prepared, "source")))
cat("  ══════════════════════════════════════════\n")

cat("[01_data_prep] Complete.\n")
