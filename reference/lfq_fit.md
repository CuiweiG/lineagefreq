# Pipe-friendly model fitting

Enables tidyverse-style chaining:
`data |> lfq_fit("mlr") |> lfq_forecast(28) |> lfq_score()`

## Usage

``` r
lfq_fit(data, engine = "mlr", ...)
```

## Arguments

- data:

  An
  [lfq_data](https://cuiweig.github.io/lineagefreq/reference/lfq_data.md)
  object.

- engine:

  Engine name. Default `"mlr"`.

- ...:

  Passed to
  [`fit_model()`](https://cuiweig.github.io/lineagefreq/reference/fit_model.md).

## Value

An `lfq_fit` object.

## Examples

``` r
data(sarscov2_us_2022)
sarscov2_us_2022 |>
  lfq_data(lineage = variant, date = date, count = count, total = total) |>
  lfq_fit("mlr") |>
  lfq_advantage(generation_time = 5)
#> # A tibble: 5 × 6
#>   lineage estimate lower upper type        pivot
#>   <chr>      <dbl> <dbl> <dbl> <chr>       <chr>
#> 1 BA.1        1     1     1    relative_Rt BA.1 
#> 2 BA.2        1.18  1.18  1.18 relative_Rt BA.1 
#> 3 BA.4/5      1.33  1.33  1.33 relative_Rt BA.1 
#> 4 BQ.1        1.29  1.28  1.29 relative_Rt BA.1 
#> 5 Other       1.11  1.11  1.11 relative_Rt BA.1 
```
