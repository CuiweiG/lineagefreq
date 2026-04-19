# Detection horizon for an emerging variant

Given current sequencing capacity and a hypothetical variant at initial
prevalence \\p_0\\ growing at rate \\\delta\\, estimates the number of
surveillance periods (weeks) until the probability of detection exceeds
a specified threshold.

## Usage

``` r
detection_horizon(
  initial_prev = 0.001,
  growth_rate = 1.3,
  n_per_period = 500L,
  n_periods = 26L,
  detection_threshold = 1L,
  confidence = 0.95
)
```

## Arguments

- initial_prev:

  Initial prevalence of the emerging variant. Default 0.001 (0.1
  percent).

- growth_rate:

  Per-week multiplicative growth rate on the frequency scale. Default
  1.3 (30 percent per week).

- n_per_period:

  Sequences collected per surveillance period. Default 500.

- n_periods:

  Maximum number of periods to evaluate. Default 26 (6 months of weekly
  data).

- detection_threshold:

  Minimum number of sequences of the target lineage required to declare
  detection. Default 1.

- confidence:

  Detection probability threshold. Default 0.95.

## Value

A tibble with columns `period`, `prevalence`, `detection_prob`
(cumulative detection probability), and `detected` (logical). An
attribute `weeks_to_detection` contains the first period where detection
probability exceeds the confidence threshold, or `NA` if not reached.

## Details

At each period \\t\\, the variant prevalence is modelled as logistic
growth: \\p(t) = p_0 \cdot r^t / (1 - p_0 + p_0 \cdot r^t)\\ where \\r\\
is the per-period growth rate. The probability of detecting at least
\\k\\ sequences in \\n\\ draws at prevalence \\p\\ is \\1 -
F\_{\text{Binom}}(k-1; n, p)\\. Cumulative detection is \\1 -
\prod\_{\tau=1}^{t}(1 - P\_\tau)\\.

## See also

[`sequencing_power`](https://cuiweig.github.io/lineagefreq/reference/sequencing_power.md)
for static sample size calculations,
[`alert_threshold`](https://cuiweig.github.io/lineagefreq/reference/alert_threshold.md)
for sequential detection.

## Examples

``` r
dh <- detection_horizon(initial_prev = 0.005, growth_rate = 1.2,
  n_per_period = 300)
dh
#> # A tibble: 26 × 4
#>    period prevalence detection_prob detected
#>     <int>      <dbl>          <dbl> <lgl>   
#>  1      1    0.00599          0.835 FALSE   
#>  2      2    0.00718          0.981 TRUE    
#>  3      3    0.00861          0.999 TRUE    
#>  4      4    0.0103           1.000 TRUE    
#>  5      5    0.0123           1.000 TRUE    
#>  6      6    0.0148           1.000 TRUE    
#>  7      7    0.0177           1.000 TRUE    
#>  8      8    0.0212           1.000 TRUE    
#>  9      9    0.0253           1     TRUE    
#> 10     10    0.0302           1     TRUE    
#> # ℹ 16 more rows
```
