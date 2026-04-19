# Information-Theoretic Surveillance Optimization

## The resource allocation problem

Genomic surveillance programmes operate under fixed budgets.
Laboratories have finite sequencing capacity, and that capacity must be
distributed across regions, time periods, and sampling strategies.
Static allocation rules (equal per region, proportional to population)
ignore the information content of each additional sequence. Optimal
allocation requires knowing where the next sequence will reduce
uncertainty the most.

`lineagefreq` v0.5.0 introduces five functions that address this problem
from an information-theoretic and decision-theoretic perspective.

## Expected Value of Information

[`surveillance_value()`](https://cuiweig.github.io/lineagefreq/reference/surveillance_value.md)
quantifies the marginal decision value of sequencing additional samples,
given the current posterior uncertainty of variant frequency estimates:

``` r
library(lineagefreq)
sim <- simulate_dynamics(n_lineages = 4, n_timepoints = 15,
  advantages = c("A" = 1.4, "B" = 1.1, "C" = 0.8), seed = 1)
fit <- fit_model(sim, engine = "mlr")

ev <- surveillance_value(fit, n_current = 500)
ev
plot(ev)
```

The EVOI curve shows diminishing returns: the first 50 additional
sequences reduce estimation variance substantially, while the next 50
contribute less. This informs budget decisions — at some point,
additional sequencing yields negligible improvement in frequency
estimates, and resources are better directed elsewhere.

## Adaptive allocation via Thompson sampling

[`adaptive_design()`](https://cuiweig.github.io/lineagefreq/reference/adaptive_design.md)
goes beyond static allocation by reallocating resources across regions
at each surveillance round, directing effort toward strata where
uncertainty reduction has the highest decision value:

``` r
ad <- adaptive_design(sim, capacity = 200, n_rounds = 10,
  strategy = "thompson", seed = 42)
ad
plot(ad)
```

Thompson sampling naturally balances exploration (sampling uncertain
regions) and exploitation (sampling where variants are most prevalent).
The alternative UCB strategy provides a deterministic upper confidence
bound approach.

This contrasts with the static Neyman allocation implemented in the
`survinger` package, which optimises for a single time point. Adaptive
allocation responds to the evolving variant landscape and is better
suited to multi-round surveillance campaigns.

## Detection horizon

[`detection_horizon()`](https://cuiweig.github.io/lineagefreq/reference/detection_horizon.md)
answers the prospective question: given current sequencing capacity and
a hypothetical new variant at initial prevalence $p_{0}$ growing at rate
$r$, how many weeks until reliable detection?

``` r
dh <- detection_horizon(initial_prev = 0.001, growth_rate = 1.3,
  n_per_period = 300, confidence = 0.95)
dh
wtd <- attr(dh, "weeks_to_detection")
```

The function accounts for multinomial sampling under logistic growth and
computes cumulative detection probability across surveillance periods.
This extends
[`sequencing_power()`](https://cuiweig.github.io/lineagefreq/reference/sequencing_power.md),
which answers the static question (how many sequences for a single
period), by adding the temporal dimension.

## Sequential detection with controlled false alarms

[`alert_threshold()`](https://cuiweig.github.io/lineagefreq/reference/alert_threshold.md)
implements sequential testing for emerging variants, controlling the
false alarm rate while minimising detection delay:

``` r
alerts <- alert_threshold(sim, method = "sprt",
  alpha = 0.05, delta_1 = 0.03)
alerts
```

Two methods are available:

- **SPRT** (Wald, 1945): Accumulates log-likelihood ratios between the
  null (stable frequency) and alternative (growing). Stops when evidence
  crosses either boundary. Controls both false alarm ($\alpha$) and
  missed detection ($\beta$) rates.

- **CUSUM**: Accumulates positive deviations from expected frequency.
  Simpler to implement but does not control missed detection probability
  explicitly.

The SPRT is particularly suited to surveillance settings where false
alarms are costly (triggering unnecessary public health responses) but
delayed detection is also costly (missing an emerging variant of
concern).

## Surveillance dashboard

[`surveillance_dashboard()`](https://cuiweig.github.io/lineagefreq/reference/surveillance_dashboard.md)
combines current landscape, growth advantages, detection power, and
calibration diagnostics into a single multi-panel display for weekly
reporting:

``` r
bt <- backtest(sim, engines = "mlr", horizons = c(7, 14),
  min_train = 42)
panels <- surveillance_dashboard(fit, sim, bt = bt)
```

The dashboard is designed for programme managers, not statisticians.
Each panel answers a specific operational question: What is circulating?
How fast is it growing? Can we detect rare variants? Are our forecasts
reliable?

## Comparison with existing approaches

No other genomic surveillance R package provides this combination of
capabilities. Static sample size calculations (phylosamp, survinger) do
not adapt to evolving dynamics. Bayesian nowcasting tools (epinowcast)
do not address resource allocation. CDC Variant Nowcast Hub produces
forecasts but not allocation recommendations or sequential detection.

The information-theoretic framework implemented here treats surveillance
as a sequential decision problem rather than a one-shot estimation
exercise.

## References

- Wald A (1945). Sequential tests of statistical hypotheses. *Annals of
  Mathematical Statistics*, 16(2), 117–186.
- Page ES (1954). Continuous inspection schemes. *Biometrika*, 41(1/2),
  100–115.
