# Test if an object is an lfq_data object

Test if an object is an lfq_data object

## Usage

``` r
is_lfq_data(x)
```

## Arguments

- x:

  Object to test.

## Value

Logical scalar.

## Examples

``` r
d <- data.frame(date = Sys.Date(), lineage = "A", count = 10)
x <- lfq_data(d, lineage = lineage, date = date, count = count)
is_lfq_data(x)
#> [1] TRUE
is_lfq_data(d)
#> [1] FALSE
```
