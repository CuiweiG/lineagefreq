# Fit model with Deep Mutational Scanning priors

A penalised multinomial logistic regression engine that incorporates
Deep Mutational Scanning (DMS) escape scores as informative priors on
variant fitness. This is valuable for early-emergence scenarios where a
new lineage has few observed sequences but laboratory-measured
phenotypic data (e.g., ACE2 binding affinity, antibody escape) are
available.

## Usage

``` r
fit_dms_prior(data, dms_scores, lambda = 1, pivot = NULL, ci_level = 0.95)
```

## Arguments

- data:

  An `lfq_data` object.

- dms_scores:

  Named numeric vector of DMS-derived fitness priors. Names correspond
  to lineage identifiers. Values are on the log growth rate scale
  (positive = fitter than average). Lineages not in the vector receive a
  prior of 0.

- lambda:

  Regularisation strength (penalty weight). Default 1.0. Larger values
  pull estimates more strongly toward the DMS prior. At `lambda = 0`,
  the result is identical to the standard MLR engine.

- pivot:

  Reference lineage name. Default `NULL` (automatic selection).

- ci_level:

  Confidence level. Default 0.95.

## Value

An `lfq_fit` object compatible with all downstream functions
(`forecast`, `growth_advantage`, etc.).

## Details

The approach uses penalised maximum likelihood where the penalty is
proportional to the squared difference between the estimated growth rate
and the DMS-derived prior. This implements an empirical Bayes shrinkage:
with abundant data, the penalty has little effect; with sparse data,
estimates are pulled toward the DMS prior.

The penalised log-likelihood is: \$\$\ell\_{\text{pen}}(\alpha, \delta)
= \ell(\alpha, \delta) - \frac{\lambda}{2} \sum_v (\delta_v -
\mu_v)^2\$\$ where \\\ell\\ is the standard multinomial log-likelihood,
\\\delta_v\\ is the growth rate for lineage \\v\\, and \\\mu_v\\ is the
DMS prior. The Hessian is adjusted accordingly, ensuring correct
confidence interval widths.

## References

Dadonaite B, Crawford KHD, Radford CE, et al. (2023). A pseudovirus
system enables deep mutational scanning of the full SARS-CoV-2 spike.
*Cell*, 186(6), 1263–1278.
[doi:10.1016/j.cell.2023.02.001](https://doi.org/10.1016/j.cell.2023.02.001)

Bloom JD, Neher RA (2023). Fitness effects of mutations to SARS-CoV-2
proteins. *Virus Evolution*, 9(2), vead055.
[doi:10.1093/ve/vead055](https://doi.org/10.1093/ve/vead055)

## See also

[`fit_model`](https://cuiweig.github.io/lineagefreq/reference/fit_model.md)
for the standard MLR engine.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.3, "B" = 0.9),
  n_timepoints = 8, total_per_tp = 100, seed = 1)
# DMS suggests lineage A has fitness advantage
dms <- c("A" = 0.04, "B" = -0.02)
fit_dms <- fit_dms_prior(sim, dms_scores = dms, lambda = 2)
growth_advantage(fit_dms)
#> # A tibble: 3 × 6
#>   lineage estimate  lower  upper type        pivot
#>   <chr>      <dbl>  <dbl>  <dbl> <chr>       <chr>
#> 1 A         0.227   0.150 0.304  growth_rate ref  
#> 2 B        -0.0882 -0.186 0.0100 growth_rate ref  
#> 3 ref       0       0     0      growth_rate ref  
# }
```
