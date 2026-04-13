# Fit a lineage frequency model

Unified interface for modeling lineage frequency dynamics. Supports
multiple engines that share the same input/output contract.

## Usage

``` r
fit_model(
  data,
  engine = c("mlr", "hier_mlr", "piantham", "fga", "garw"),
  pivot = NULL,
  horizon = 28L,
  ci_level = 0.95,
  ...
)
```

## Arguments

- data:

  An
  [lfq_data](https://CuiweiG.github.io/lineagefreq/reference/lfq_data.md)
  object.

- engine:

  Model engine to use:

  - `"mlr"` (default): Multinomial logistic regression. Fast,
    frequentist, no external dependencies.

  - `"hier_mlr"`: Hierarchical MLR with partial pooling across
    locations. Requires `.location` column in data.

  - `"piantham"`: Piantham approximation converting MLR growth rates to
    relative reproduction numbers. Requires `generation_time` argument.

  - `"fga"`: Fixed growth advantage model (Bayesian via CmdStan).
    Requires 'CmdStan'; check with
    [`lfq_stan_available()`](https://CuiweiG.github.io/lineagefreq/reference/lfq_stan_available.md).

  - `"garw"`: Growth advantage random walk model (Bayesian via CmdStan).
    Allows fitness to change over time.

- pivot:

  Reference lineage name. Growth rates are reported relative to this
  lineage (fixed at 0). Default: the lineage with the highest count at
  the earliest time point.

- horizon:

  Forecast horizon in days (stored for later use by
  [`forecast()`](https://CuiweiG.github.io/lineagefreq/reference/forecast.md)).
  Default 28.

- ci_level:

  Confidence level for intervals. Default 0.95.

- ...:

  Engine-specific arguments passed to the internal engine function. For
  `engine = "mlr"`: `window`, `ci_method`, `laplace_smooth`. For
  `engine = "piantham"`: `generation_time` (required). For
  `engine = "hier_mlr"`: `shrinkage_method`.

## Value

An `lfq_fit` object (S3 class), a list containing:

- engine:

  Engine name (character).

- growth_rates:

  Named numeric vector of growth rates per `time_scale` days (pivot =
  0).

- intercepts:

  Named numeric vector of intercepts.

- pivot:

  Name of pivot lineage.

- lineages:

  Character vector of all lineage names.

- fitted_values:

  Tibble of fitted frequencies.

- residuals:

  Tibble with observed, fitted, Pearson residuals.

- vcov_matrix:

  Variance-covariance matrix.

- loglik, aic, bic:

  Model fit statistics.

- nobs, n_timepoints, df:

  Sample and model sizes.

- ci_level, horizon:

  As specified.

- call:

  The matched call.

## Details

The MLR engine models the frequency of lineage \\v\\ at time \\t\\ as:
\$\$p_v(t) = \frac{\exp(\alpha_v + \delta_v t)}{\sum_k \exp(\alpha_k +
\delta_k t)}\$\$ where \\\alpha_v\\ is the intercept, \\\delta_v\\ is
the growth rate per `time_scale` days (default 7), and the pivot lineage
has \\\alpha = \delta = 0\\. Parameters are estimated by maximum
likelihood via `optim(method = "BFGS")` with the Hessian used for
asymptotic Wald confidence intervals.

The constant growth rate assumption is appropriate for monotonic variant
replacement periods (typically 2–4 months). For longer periods or
non-monotonic dynamics, use the `window` argument to restrict the
estimation window, or consider the `"garw"` engine which allows
time-varying growth advantages.

## References

Abousamra E, Figgins M, Bedford T (2024). Fitness models provide
accurate short-term forecasts of SARS-CoV-2 variant frequency. *PLoS
Computational Biology*, 20(9):e1012443.
[doi:10.1371/journal.pcbi.1012443](https://doi.org/10.1371/journal.pcbi.1012443)

Piantham C, Linton NM, Nishiura H (2022). Predicting the trajectory of
replacements of SARS-CoV-2 variants using relative reproduction numbers.
*Viruses*, 14(11):2556.
[doi:10.3390/v14112556](https://doi.org/10.3390/v14112556)

## See also

[`growth_advantage()`](https://CuiweiG.github.io/lineagefreq/reference/growth_advantage.md)
to extract fitness estimates,
[`forecast()`](https://CuiweiG.github.io/lineagefreq/reference/forecast.md)
for frequency prediction,
[`backtest()`](https://CuiweiG.github.io/lineagefreq/reference/backtest.md)
for rolling-origin evaluation.

## Examples

``` r
sim <- simulate_dynamics(
  n_lineages = 3,
  advantages = c("JN.1" = 1.3, "KP.3" = 0.9),
  n_timepoints = 15, seed = 42
)
fit <- fit_model(sim, engine = "mlr")
fit
#> Lineage frequency model (mlr)
#> 3 lineages, 15 time points
#> Date range: 2026-04-13 to 2026-07-20
#> Pivot: "KP.3"
#> 
#> Growth rates (per 7-day unit):
#>   ↑ JN.1: 0.3638
#>   ↑ ref: 0.1122
#> 
#> AIC: 9450; BIC: 9453
```
