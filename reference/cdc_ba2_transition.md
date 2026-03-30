# CDC SARS-CoV-2 variant proportions: BA.1 to BA.2 transition (US, 2022)

Real surveillance data covering the Omicron BA.1 to BA.2 variant
replacement in the United States, December 2021 through June 2022. This
is one of the best-documented variant replacement events and serves as
an independent validation dataset.

## Usage

``` r
cdc_ba2_transition
```

## Format

A data frame with 150 rows and 4 columns:

- date:

  Biweek ending date (Date).

- lineage:

  Lineage name: BA.1, BA.2, BA.2.12.1, BA.4/5, Other.

- count:

  Approximate sequence count per biweek (integer).

- proportion:

  CDC weighted proportion estimate (numeric).

## Source

CDC COVID Data Tracker (data.cdc.gov, public domain).

## References

Lyngse FP, et al. (2022). Household transmission of SARS-CoV-2 Omicron
variant of concern subvariants BA.1 and BA.2 in Denmark. *Nature
Communications*, 13:5760.
[doi:10.1038/s41467-022-33498-0](https://doi.org/10.1038/s41467-022-33498-0)

## Examples

``` r
data(cdc_ba2_transition)
vd <- lfq_data(cdc_ba2_transition,
               date = date, lineage = lineage, count = count)
fit <- fit_model(vd, engine = "mlr", pivot = "BA.1")
# BA.2 Rt ~ 1.34 (consistent with published estimates)
growth_advantage(fit, type = "relative_Rt", generation_time = 3.2)
#> # A tibble: 5 × 6
#>   lineage   estimate lower upper type        pivot
#>   <chr>        <dbl> <dbl> <dbl> <chr>       <chr>
#> 1 BA.1         1     1     1     relative_Rt BA.1 
#> 2 BA.2         1.34  1.34  1.34  relative_Rt BA.1 
#> 3 BA.2.12.1    1.56  1.55  1.56  relative_Rt BA.1 
#> 4 BA.4/5       2.01  2.00  2.02  relative_Rt BA.1 
#> 5 Other        0.707 0.704 0.711 relative_Rt BA.1 
```
