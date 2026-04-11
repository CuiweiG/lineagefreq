# Convert lfq_data to long-format tibble

Convert lfq_data to long-format tibble

## Usage

``` r
# S3 method for class 'lfq_data'
as.data.frame(x, ...)
```

## Arguments

- x:

  An
  [lfq_data](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
  object.

- ...:

  Ignored.

## Value

A tibble with all columns.

## Examples

``` r
data(sarscov2_us_2022)
x <- lfq_data(sarscov2_us_2022, lineage = variant,
              date = date, count = count, total = total)
as.data.frame(x)
#> # A tibble: 200 × 7
#>    .date      .lineage .count total .total    .freq .reliable
#>    <date>     <chr>     <int> <int>  <int>    <dbl> <lgl>    
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
