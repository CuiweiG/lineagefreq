# Coerce to lfq_data

Generic function to convert various data formats into
[lfq_data](https://cuiweig.github.io/lineagefreq/reference/lfq_data.md)
objects. Methods can be defined for specific input classes to enable
seamless interoperability with other genomic surveillance packages.

## Usage

``` r
as_lfq_data(x, ...)

# S3 method for class 'lfq_data'
as_lfq_data(x, ...)

# S3 method for class 'data.frame'
as_lfq_data(x, ...)
```

## Arguments

- x:

  An object to coerce.

- ...:

  Additional arguments passed to methods. For the `data.frame` method,
  these are passed to
  [`lfq_data()`](https://cuiweig.github.io/lineagefreq/reference/lfq_data.md).

## Value

An
[lfq_data](https://cuiweig.github.io/lineagefreq/reference/lfq_data.md)
object.

An
[lfq_data](https://cuiweig.github.io/lineagefreq/reference/lfq_data.md)
object.

An
[lfq_data](https://cuiweig.github.io/lineagefreq/reference/lfq_data.md)
object.

## See also

[`lfq_data()`](https://cuiweig.github.io/lineagefreq/reference/lfq_data.md)
for the primary constructor.

## Examples

``` r
df <- data.frame(
  date    = rep(as.Date("2024-01-01") + c(0, 7), each = 2),
  lineage = rep(c("A", "B"), 2),
  count   = c(80, 20, 60, 40)
)
x <- as_lfq_data(df, lineage = lineage, date = date, count = count)
x
#> 
#> ── Lineage frequency data 
#> 2 lineages, 2 time points
#> Date range: 2024-01-01 to 2024-01-08
#> Lineages: "A, B"
#> 
#> # A tibble: 4 × 6
#>   .date      .lineage .count .total .freq .reliable
#> * <date>     <chr>     <int>  <int> <dbl> <lgl>    
#> 1 2024-01-01 A            80    100   0.8 TRUE     
#> 2 2024-01-01 B            20    100   0.2 TRUE     
#> 3 2024-01-08 A            60    100   0.6 TRUE     
#> 4 2024-01-08 B            40    100   0.4 TRUE     
```
