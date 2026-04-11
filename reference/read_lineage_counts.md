# Read lineage count data from a CSV file

A convenience wrapper for reading surveillance count data from CSV files
into
[lfq_data](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
format. Expects a file with at least columns for date, lineage, and
count.

## Usage

``` r
read_lineage_counts(
  file,
  date = "date",
  lineage = "lineage",
  count = "count",
  ...
)
```

## Arguments

- file:

  Path to CSV file.

- date:

  Name of the date column. Default `"date"`.

- lineage:

  Name of the lineage column. Default `"lineage"`.

- count:

  Name of the count column. Default `"count"`.

- ...:

  Additional arguments passed to
  [`lfq_data()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md).

## Value

An
[lfq_data](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
object.

## Examples

``` r
# Read the bundled example CSV
f <- system.file("extdata", "example_counts.csv",
                 package = "lineagefreq")
x <- read_lineage_counts(f)
x
#> 
#> ── Lineage frequency data 
#> 5 lineages, 3 time points
#> Date range: 2021-12-04 to 2021-12-18
#> Lineages: "BA.1, BA.2, BA.2.12.1, BA.4/5, Other"
#> 
#> # A tibble: 15 × 6
#>    .date      .lineage  .count .total    .freq .reliable
#>  * <date>     <chr>      <int>  <int>    <dbl> <lgl>    
#>  1 2021-12-04 BA.1          61   8000 0.00762  TRUE     
#>  2 2021-12-04 BA.2           0   8000 0        TRUE     
#>  3 2021-12-04 BA.2.12.1      0   8000 0        TRUE     
#>  4 2021-12-04 BA.4/5         0   8000 0        TRUE     
#>  5 2021-12-04 Other       7939   8000 0.992    TRUE     
#>  6 2021-12-11 BA.1         529   8000 0.0661   TRUE     
#>  7 2021-12-11 BA.2           1   8000 0.000125 TRUE     
#>  8 2021-12-11 BA.2.12.1      0   8000 0        TRUE     
#>  9 2021-12-11 BA.4/5         0   8000 0        TRUE     
#> 10 2021-12-11 Other       7470   8000 0.934    TRUE     
#> 11 2021-12-18 BA.1        3022   8000 0.378    TRUE     
#> 12 2021-12-18 BA.2           2   8000 0.00025  TRUE     
#> 13 2021-12-18 BA.2.12.1      0   8000 0        TRUE     
#> 14 2021-12-18 BA.4/5         0   8000 0        TRUE     
#> 15 2021-12-18 Other       4976   8000 0.622    TRUE     
```
