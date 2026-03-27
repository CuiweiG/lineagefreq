# Create a lineage frequency data object

Validates, structures, and annotates lineage count data for downstream
modeling and analysis. This is the entry point for all lineagefreq
workflows.

## Usage

``` r
lfq_data(
  data,
  lineage,
  date,
  count,
  total = NULL,
  location = NULL,
  min_total = 10L
)
```

## Arguments

- data:

  A data frame containing at minimum columns for lineage identity, date,
  and count.

- lineage:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Column containing lineage/variant identifiers (character or factor).

- date:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Column containing collection dates (`Date` class or parseable
  character).

- count:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Column containing sequence counts (non-negative integers).

- total:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optional column of total sequences per date-location. If `NULL`,
  computed as the sum of `count` per group.

- location:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optional column for geographic stratification.

- min_total:

  Minimum total count per time point. Time points below this are flagged
  as unreliable. Default 10.

## Value

An `lfq_data` object (a tibble subclass) with standardized columns:

- `.lineage`:

  Lineage identifier (character).

- `.date`:

  Collection date (Date).

- `.count`:

  Sequence count (integer).

- `.total`:

  Total sequences at this time point (integer).

- `.freq`:

  Observed frequency (numeric).

- `.reliable`:

  Logical; `TRUE` if `.total >= min_total`.

- `.location`:

  Location, if provided (character).

All original columns are preserved.

## Details

Performs the following validation and processing:

1.  Checks that all required columns exist and have correct types.

2.  Coerces character dates to Date via ISO 8601 parsing.

3.  Ensures counts are non-negative integers.

4.  Replaces NA counts with 0 (with warning).

5.  Aggregates duplicate lineage-date rows by summing (with warning).

6.  Computes per-time-point totals and frequencies.

7.  Flags time points below `min_total` as unreliable.

8.  Sorts by date ascending, then lineage alphabetically.

## Examples

``` r
d <- data.frame(
  date = rep(seq(as.Date("2024-01-01"), by = "week",
                 length.out = 8), each = 3),
  lineage = rep(c("JN.1", "KP.3", "Other"), 8),
  n = c(5, 2, 93, 12, 5, 83, 28, 11, 61, 50, 20, 30,
        68, 18, 14, 80, 12, 8, 88, 8, 4, 92, 5, 3)
)
x <- lfq_data(d, lineage = lineage, date = date, count = n)
x
#> 
#> ── Lineage frequency data 
#> 3 lineages, 8 time points
#> Date range: 2024-01-01 to 2024-02-19
#> Lineages: "JN.1, KP.3, Other"
#> 
#> # A tibble: 24 × 6
#>    .date      .lineage .count .total .freq .reliable
#>  * <date>     <chr>     <int>  <int> <dbl> <lgl>    
#>  1 2024-01-01 JN.1          5    100  0.05 TRUE     
#>  2 2024-01-01 KP.3          2    100  0.02 TRUE     
#>  3 2024-01-01 Other        93    100  0.93 TRUE     
#>  4 2024-01-08 JN.1         12    100  0.12 TRUE     
#>  5 2024-01-08 KP.3          5    100  0.05 TRUE     
#>  6 2024-01-08 Other        83    100  0.83 TRUE     
#>  7 2024-01-15 JN.1         28    100  0.28 TRUE     
#>  8 2024-01-15 KP.3         11    100  0.11 TRUE     
#>  9 2024-01-15 Other        61    100  0.61 TRUE     
#> 10 2024-01-22 JN.1         50    100  0.5  TRUE     
#> # ℹ 14 more rows
```
