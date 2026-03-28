# data-raw/

This directory contains scripts that reproduce the package's built-in
datasets from the original CDC source data.

## Data provenance

| Dataset | Source | License | Script |
|---------|--------|---------|--------|
| `cdc_sarscov2_jn1` | CDC Variant Proportions (jr58-6ysp) | Public domain (17 USC 105) | `prepare_cdc_data.R` |
| `cdc_ba2_transition` | CDC Variant Proportions (jr58-6ysp) | Public domain (17 USC 105) | `prepare_ba2_data.R` |
| `sarscov2_us_2022` | Simulated | — | `prepare_datasets.R` |
| `influenza_h3n2` | Simulated | — | `prepare_datasets.R` |

## Reproducing from source

1. Download the full CDC CSV (≈800 MB):
   ```
   https://data.cdc.gov/api/views/jr58-6ysp/rows.csv?accessType=DOWNLOAD
   ```

2. Place as `SARS-CoV-2_Variant_Proportions.csv` in this directory.

3. Run the scripts:
   ```r
   source("data-raw/prepare_cdc_data.R")
   source("data-raw/prepare_ba2_data.R")
   ```

## Processing steps

Each script:
1. Reads the raw CDC CSV (4.4M rows)
2. Filters to national-level weighted estimates for the target period
3. Aggregates variant categories (collapsing sublineages)
4. Converts CDC's proportion estimates to approximate integer counts
   (preserving multinomial sampling properties)
5. Saves as `.rda` with `compress = "xz"`

The `count` column represents approximate sequence counts reconstructed
from CDC's published proportions. The `proportion` column retains the
original CDC estimates for reference.

## References

Ma KC, et al. (2024). Genomic Surveillance for SARS-CoV-2 Variants.
*MMWR*, 73(42):941–948. doi:10.15585/mmwr.mm7342a1
