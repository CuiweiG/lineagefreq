# Simulated SARS-CoV-2 variant frequency data (US, 2022)

A simulated dataset of weekly SARS-CoV-2 variant sequence counts for the
United States in 2022. Includes the BA.1 to BA.2 to BA.4/5 to BQ.1
transition dynamics. This is simulated data (not real GISAID data) to
avoid license restrictions while preserving realistic statistical
properties.

## Usage

``` r
sarscov2_us_2022
```

## Format

A data frame with 200 rows and 4 columns:

- date:

  Collection date (Date, weekly).

- variant:

  Variant name (character): BA.1, BA.2, BA.4/5, BQ.1, Other.

- count:

  Number of sequences assigned to this variant in this week (integer).

- total:

  Total sequences for this week (integer).

## Source

Simulated based on parameters from published CDC MMWR genomic
surveillance reports and 'Nextstrain' public data.

## Examples

``` r
data(sarscov2_us_2022)
x <- lfq_data(sarscov2_us_2022, lineage = variant,
              date = date, count = count, total = total)
x
#> 
#> ── Lineage frequency data 
#> 5 lineages, 40 time points
#> Date range: 2022-01-08 to 2022-10-08
#> Lineages: "BA.1, BA.2, BA.4/5, BQ.1, Other"
#> 
#> # A tibble: 200 × 7
#>    .date      .lineage .count total .total    .freq .reliable
#>  * <date>     <chr>     <int> <int>  <int>    <dbl> <lgl>    
#>  1 2022-01-08 BA.1      11802 14929  14929 0.791    TRUE     
#>  2 2022-01-08 BA.2        565 14929  14929 0.0378   TRUE     
#>  3 2022-01-08 BA.4/5        4 14929  14929 0.000268 TRUE     
#>  4 2022-01-08 BQ.1          0 14929  14929 0        TRUE     
#>  5 2022-01-08 Other      2558 14929  14929 0.171    TRUE     
#>  6 2022-01-15 BA.1      10342 13672  13672 0.756    TRUE     
#>  7 2022-01-15 BA.2        608 13672  13672 0.0445   TRUE     
#>  8 2022-01-15 BA.4/5        3 13672  13672 0.000219 TRUE     
#>  9 2022-01-15 BQ.1          0 13672  13672 0        TRUE     
#> 10 2022-01-15 Other      2719 13672  13672 0.199    TRUE     
#> # ℹ 190 more rows
```
