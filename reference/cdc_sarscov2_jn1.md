# CDC SARS-CoV-2 variant proportions: JN.1 emergence (US, 2023-2024)

Real surveillance data from the CDC's national genomic surveillance
program covering the emergence and dominance of the SARS-CoV-2 JN.1
lineage in the United States, October 2023 through June 2024.

## Usage

``` r
cdc_sarscov2_jn1
```

## Format

A data frame with 171 rows and 4 columns:

- date:

  Biweek ending date (Date).

- lineage:

  Lineage name (character): JN.1, XBB.1.5, EG.5.1, HV.1, HK.3, BA.2.86,
  KP.2, KP.3, JN.1.11.1, Other.

- count:

  Approximate sequence count per biweek (integer).

- proportion:

  CDC weighted proportion estimate (numeric).

## Source

CDC COVID Data Tracker, SARS-CoV-2 Variant Proportions. Dataset ID:
jr58-6ysp.
<https://data.cdc.gov/Laboratory-Surveillance/SARS-CoV-2-Variant-Proportions/jr58-6ysp>

Public domain (U.S. Government Work, 17 USC 105).

## Details

Derived from CDC's published weighted variant proportion estimates.
Approximate biweekly sequence counts were reconstructed from proportions
using a reference total of 8,000 sequences per period. The original
proportions are retained in the `proportion` column.

## References

Ma KC, et al. (2024). Genomic Surveillance for SARS-CoV-2 Variants.
*MMWR*, 73(42):941–948.
[doi:10.15585/mmwr.mm7342a1](https://doi.org/10.15585/mmwr.mm7342a1)

## Examples

``` r
data(cdc_sarscov2_jn1)
vd <- lfq_data(cdc_sarscov2_jn1,
               date = date, lineage = lineage, count = count)
fit <- fit_model(vd, engine = "mlr")
growth_advantage(fit, type = "relative_Rt", generation_time = 5)
#> # A tibble: 9 × 6
#>   lineage   estimate lower upper type        pivot
#>   <chr>        <dbl> <dbl> <dbl> <chr>       <chr>
#> 1 BA.2.86      0.961 0.957 0.965 relative_Rt Other
#> 2 HK.3         0.900 0.896 0.904 relative_Rt Other
#> 3 HV.1         0.898 0.896 0.900 relative_Rt Other
#> 4 JN.1         1.01  1.01  1.01  relative_Rt Other
#> 5 JN.1.11.1    1.11  1.10  1.12  relative_Rt Other
#> 6 KP.2         1.16  1.15  1.16  relative_Rt Other
#> 7 KP.3         1.26  1.25  1.27  relative_Rt Other
#> 8 Other        1     1     1     relative_Rt Other
#> 9 XBB.1.5      0.860 0.850 0.870 relative_Rt Other
```
