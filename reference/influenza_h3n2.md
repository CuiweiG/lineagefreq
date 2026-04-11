# Simulated influenza A/H3N2 clade frequency data

A simulated dataset of weekly influenza A/H3N2 clade sequence counts
over a single Northern Hemisphere season (24 weeks).

## Usage

``` r
influenza_h3n2
```

## Format

A data frame with 96 rows and 4 columns:

- date:

  Collection date (Date, weekly).

- clade:

  Clade name (character).

- count:

  Number of sequences (integer).

- total:

  Total sequences for this week (integer).

## Source

Simulated based on 'Nextstrain' influenza clade dynamics.

## Examples

``` r
data(influenza_h3n2)
x <- lfq_data(influenza_h3n2, lineage = clade,
              date = date, count = count, total = total)
x
#> 
#> ── Lineage frequency data 
#> 4 lineages, 24 time points
#> Date range: 2023-10-01 to 2024-03-10
#> Lineages: "2a1b.1, 2a1b.2a.1, 2a1b.2a.2, Other"
#> 
#> # A tibble: 96 × 7
#>    .date      .lineage  .count total .total  .freq .reliable
#>  * <date>     <chr>      <int> <int>  <int>  <dbl> <lgl>    
#>  1 2023-10-01 2a1b.1       137   515    515 0.266  TRUE     
#>  2 2023-10-01 2a1b.2a.1    234   515    515 0.454  TRUE     
#>  3 2023-10-01 2a1b.2a.2     51   515    515 0.0990 TRUE     
#>  4 2023-10-01 Other         93   515    515 0.181  TRUE     
#>  5 2023-10-08 2a1b.1        84   303    303 0.277  TRUE     
#>  6 2023-10-08 2a1b.2a.1    123   303    303 0.406  TRUE     
#>  7 2023-10-08 2a1b.2a.2     45   303    303 0.149  TRUE     
#>  8 2023-10-08 Other         51   303    303 0.168  TRUE     
#>  9 2023-10-15 2a1b.1        95   451    451 0.211  TRUE     
#> 10 2023-10-15 2a1b.2a.1    196   451    451 0.435  TRUE     
#> # ℹ 86 more rows
```
