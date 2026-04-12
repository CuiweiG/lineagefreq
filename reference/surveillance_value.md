# Expected Value of Information for genomic surveillance

Quantifies the marginal decision value of sequencing additional samples
from each region or stratum, based on the current posterior uncertainty
of variant frequency estimates. The EVOI captures how much the expected
estimation error decreases per additional sequence, enabling
cost-effective resource allocation.

## Usage

``` r
surveillance_value(
  fit,
  n_current = 100L,
  n_additional = seq(10L, 200L, by = 10L),
  target_lineage = NULL
)
```

## Arguments

- fit:

  An `lfq_fit` object from
  [`fit_model`](https://CuiweiG.github.io/lineagefreq/reference/fit_model.md).

- n_current:

  Integer vector of current sample sizes per stratum. If a single
  integer, assumed equal across all lineages.

- n_additional:

  Integer vector of candidate additional sample sizes to evaluate.
  Default `seq(10, 200, by = 10)`.

- target_lineage:

  Optional character; compute EVOI specifically for this lineage.
  Default `NULL` evaluates across all non-pivot lineages.

## Value

An `evoi` S3 class with components:

- values:

  Tibble with `n_additional`, `evoi`, `marginal_evoi` columns.

- current_uncertainty:

  Numeric; current estimation variance (sum of growth rate SEs squared).

- target_lineage:

  Character or NULL.

## Details

The EVOI is computed as the expected reduction in mean squared
estimation error for lineage frequency when \\n\\ additional sequences
are observed. Under the multinomial likelihood with Gaussian
approximation to the posterior, the variance of the frequency estimate
scales as \\p(1-p) / n\\, and additional samples reduce this by the
factor \\n_0 / (n_0 + n\_{\text{add}})\\.

The marginal EVOI (the value of one additional sequence) is the
derivative of the EVOI curve. It decreases monotonically with sample
size, exhibiting the diminishing returns characteristic of
information-theoretic quantities.

## See also

[`adaptive_design`](https://CuiweiG.github.io/lineagefreq/reference/adaptive_design.md)
for dynamic allocation,
[`sequencing_power`](https://CuiweiG.github.io/lineagefreq/reference/sequencing_power.md)
for static sample size planning.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.3, "B" = 0.9),
  n_timepoints = 15, seed = 1)
fit <- fit_model(sim, engine = "mlr")
ev <- surveillance_value(fit, n_current = 500)
ev
#> 
#> ── Expected Value of Information 
#> Current estimation variance: 0.00032
#> 50% variance reduction at ~80 additional sequences
# }
```
