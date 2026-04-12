# Sequential detection of emerging variants

Applies a sequential probability ratio test (SPRT) or CUSUM procedure to
lineage frequency data, determining when accumulated evidence is
sufficient to declare a variant "emerging" rather than sampling noise.
Controls the false alarm rate while minimising detection delay.

## Usage

``` r
alert_threshold(
  data,
  method = c("sprt", "cusum"),
  alpha = 0.05,
  beta = 0.1,
  delta_0 = 0,
  delta_1 = 0.03,
  threshold = 5
)
```

## Arguments

- data:

  An `lfq_data` object.

- method:

  Detection method: `"sprt"` (default) for the sequential probability
  ratio test, or `"cusum"` for the cumulative sum control chart.

- alpha:

  False alarm probability. Default 0.05.

- beta:

  Missed detection probability. Default 0.10.

- delta_0:

  Null hypothesis growth rate (no emergence). Default 0 (frequency is
  stable).

- delta_1:

  Alternative hypothesis growth rate (emergence). Default 0.03 (3
  percent per-week increase on logit scale).

- threshold:

  CUSUM decision threshold. Default 5.0. Only used when
  `method = "cusum"`.

## Value

A tibble with columns `lineage`, `date`, `statistic` (log-likelihood
ratio or CUSUM value), `alert` (logical), `direction`
(emerging/declining/ stable), and `confidence` (1 - alpha).

## Details

**SPRT** (Wald, 1945) computes the log-likelihood ratio between the
alternative (lineage is growing at rate \\\delta_1\\) and the null
(frequency is stable). The test stops when the cumulative log-ratio
crosses the upper boundary \\B = \log((1-\beta)/\alpha)\\ (declare
emerging) or the lower boundary \\A = \log(\beta/(1-\alpha))\\ (declare
stable).

**CUSUM** accumulates deviations from expected frequency under the null:
\\S_t = \max(0, S\_{t-1} + (x_t - k))\\ where \\x_t\\ is the observed
frequency change and \\k\\ is the allowance (half the shift to detect).
An alert is raised when \\S_t \> h\\.

## References

Wald A (1945). Sequential tests of statistical hypotheses. *Annals of
Mathematical Statistics*, 16(2), 117–186.

## See also

[`summarize_emerging`](https://CuiweiG.github.io/lineagefreq/reference/summarize_emerging.md)
for non-sequential trend tests,
[`detection_horizon`](https://CuiweiG.github.io/lineagefreq/reference/detection_horizon.md)
for prospective power analysis.

## Examples

``` r
# \donttest{
sim <- simulate_dynamics(n_lineages = 3,
  advantages = c("A" = 1.5, "B" = 0.8),
  n_timepoints = 15, seed = 1)
alerts <- alert_threshold(sim)
alerts
#> # A tibble: 3 × 6
#>   lineage date       statistic alert direction    confidence
#>   <chr>   <date>         <dbl> <lgl> <chr>             <dbl>
#> 1 A       2026-07-19     0.222 FALSE inconclusive       0.95
#> 2 B       2026-07-19    -0.193 FALSE inconclusive       0.95
#> 3 ref     2026-07-19    -0.193 FALSE inconclusive       0.95
# }
```
