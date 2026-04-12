# Construct a population immunity landscape

Assembles time-varying population immunity estimates against each
circulating lineage from seroprevalence surveys, vaccination records, or
infection history data. The resulting object serves as input to
[`fitness_decomposition`](https://CuiweiG.github.io/lineagefreq/reference/fitness_decomposition.md)
for disentangling intrinsic transmissibility from immune escape.

## Usage

``` r
immune_landscape(
  data,
  date,
  lineage,
  immunity,
  type = "combined",
  cross_immunity = NULL
)
```

## Arguments

- data:

  A data frame with columns for date, lineage, and immunity level.

- date:

  Column name (unquoted or string) containing dates.

- lineage:

  Column name containing lineage/variant identifiers.

- immunity:

  Column name containing population-level immunity estimates
  (proportion, 0–1 scale) against each lineage.

- type:

  Character vector specifying immunity source(s): `"infection"`,
  `"vaccine"`, `"hybrid"`, or `"combined"` (default). If the data
  contain a `type` column, per-source estimates are preserved; otherwise
  all estimates are treated as combined.

- cross_immunity:

  Optional numeric matrix of cross-immunity between lineages. Rows and
  columns correspond to lineages; entry \\(i,j)\\ is the degree to which
  immunity against lineage \\i\\ protects against lineage \\j\\ (0 to
  1). Default `NULL` assumes no cross-protection data.

## Value

An `immune_landscape` object (S3 class) with components:

- estimates:

  Tibble with columns `date`, `lineage`, `immunity`, and optionally
  `type`.

- lineages:

  Character vector of lineage names.

- date_range:

  Two-element Date vector.

- cross_immunity:

  Matrix or NULL.

## Details

Immunity estimates may come from multiple sources. Seroprevalence
surveys provide direct measurements but are infrequent and
geographically sparse. Vaccination coverage data are more widely
available but do not capture waning or variant-specific escape.
Model-based reconstructions (e.g., from case and death data) can fill
gaps but introduce model dependence.

`immune_landscape()` accepts any of these as input. The critical
requirement is that immunity is expressed on a 0–1 scale representing
the proportion of the population with neutralising protection against
each lineage at each time point.

## See also

[`fitness_decomposition`](https://CuiweiG.github.io/lineagefreq/reference/fitness_decomposition.md)
for downstream analysis.

## Examples

``` r
# Simulated immunity data
imm_data <- data.frame(
  date = rep(seq(as.Date("2024-01-01"), by = "week",
                 length.out = 10), each = 3),
  lineage = rep(c("BA.5", "XBB.1.5", "JN.1"), 10),
  immunity = c(
    rep(c(0.6, 0.4, 0.1), 5),
    rep(c(0.55, 0.45, 0.25), 5))
)
il <- immune_landscape(imm_data, date = date,
  lineage = lineage, immunity = immunity)
il
#> 
#> ── Population immunity landscape 
#> 3 lineages: "BA.5, JN.1, XBB.1.5"
#> Date range: 2024-01-01 to 2024-03-04
#> Immunity types: "combined"
```
