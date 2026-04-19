# Decompose variant fitness into transmissibility and immune escape

Partitions the observed growth advantage of each lineage into two
components: intrinsic transmissibility (fitness in a fully naive
population) and immune escape (additional fitness gained from evading
population immunity). This decomposition follows the framework of
Figgins and Bedford (2025), where the effective fitness of lineage \\v\\
is modelled as \\f_v(t) = \beta_v \cdot (1 - \pi_v(t))\\ with
\\\beta_v\\ representing intrinsic transmissibility and \\\pi_v(t)\\ the
proportion of the population with neutralising immunity against \\v\\.

## Usage

``` r
fitness_decomposition(fit, landscape, generation_time)
```

## Arguments

- fit:

  An `lfq_fit` object from
  [`fit_model`](https://cuiweig.github.io/lineagefreq/reference/fit_model.md).

- landscape:

  An `immune_landscape` object from
  [`immune_landscape`](https://cuiweig.github.io/lineagefreq/reference/immune_landscape.md).

- generation_time:

  Generation time in days (required for converting growth rates to
  reproduction numbers).

## Value

A `fitness_decomposition` object (S3 class) with components:

- decomposition:

  Tibble with columns: `lineage`, `observed_advantage` (total growth
  rate), `beta` (intrinsic transmissibility), `escape_contribution`
  (immune escape component), `transmissibility_fraction` (proportion due
  to intrinsic fitness), `escape_fraction` (proportion due to immune
  escape).

- fit:

  The input lfq_fit object.

- landscape:

  The input immune_landscape.

- generation_time:

  Generation time used.

## Details

The decomposition proceeds as follows. For each non-pivot lineage, the
observed growth rate \\\delta_v\\ from the MLR fit is taken as the total
fitness difference relative to the reference. The immune escape
component is estimated from the immunity differential:
\$\$\text{escape}\_v = -\log(1 - \bar{\pi}\_v) + \log(1 -
\bar{\pi}\_{\text{ref}})\$\$ where \\\bar{\pi}\_v\\ is the mean
population immunity against lineage \\v\\ over the fitted period. The
intrinsic transmissibility component is the residual: \$\$\beta_v =
\delta_v - \text{escape}\_v\$\$

When cross-immunity data are available in the landscape object,
effective immunity against each lineage is computed as the weighted sum
of immunity sources, discounted by the cross-immunity matrix.

## References

Figgins MD, Bedford T (2025). Jointly modeling variant frequencies and
case counts to estimate relative variant severity. *medRxiv*.
[doi:10.1101/2024.12.02.24318334](https://doi.org/10.1101/2024.12.02.24318334)

## See also

[`immune_landscape`](https://cuiweig.github.io/lineagefreq/reference/immune_landscape.md)
for constructing the immunity input,
[`selective_pressure`](https://cuiweig.github.io/lineagefreq/reference/selective_pressure.md)
for population-level early warning signals.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.3, "B" = 0.9),
  n_timepoints = 15, seed = 1)
fit <- fit_model(sim, engine = "mlr")

imm_data <- data.frame(
  date = rep(unique(sim$.date), each = 3),
  lineage = rep(c("A", "B", "ref"), length(unique(sim$.date))),
  immunity = c(rep(c(0.2, 0.5, 0.4),
    length(unique(sim$.date))))
)
il <- immune_landscape(imm_data, date, lineage, immunity)
fd <- fitness_decomposition(fit, il, generation_time = 5)
fd
#> 
#> ── Fitness decomposition 
#> Pivot: ref | Generation time: 5 days
#> 
#> A: 69% intrinsic / 31% escape (delta = 0.2579)
#> B: 64% intrinsic / 36% escape (delta = -0.1061)
# }
```
