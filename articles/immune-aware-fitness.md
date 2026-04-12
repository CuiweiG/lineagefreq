# Immune-Aware Fitness Estimation

## The confounding of transmissibility and immune escape

Standard multinomial logistic regression estimates a single growth rate
per lineage. But the observed growth advantage of a new variant
conflates two distinct mechanisms: intrinsic transmissibility (the
variant’s ability to spread in a naive population) and immune escape
(its ability to evade existing population immunity). A variant with
moderate intrinsic transmissibility but strong immune escape can have
the same observed growth rate as one with high intrinsic
transmissibility and no escape at all. The public health implications of
these two scenarios are quite different.

`lineagefreq` v0.4.0 introduces tools for disentangling these
components, incorporating external immunological data, and computing
early warning signals from genomic data alone.

## Constructing the immunity landscape

The first step is encoding what is known about population-level immunity
against each circulating lineage. This may come from seroprevalence
surveys, vaccination coverage data, or model-based reconstructions.

``` r
library(lineagefreq)

# Immunity data: proportion of population with neutralising
# protection against each lineage at each time point
imm_data <- data.frame(
  date = rep(seq(as.Date("2024-01-01"), by = "week",
                 length.out = 15), each = 3),
  lineage = rep(c("JN.1", "XBB.1.5", "BA.5"), 15),
  immunity = c(
    rep(c(0.10, 0.55, 0.65), 8),
    rep(c(0.25, 0.50, 0.60), 7)
  )
)

il <- immune_landscape(imm_data, date = date,
                       lineage = lineage, immunity = immunity)
plot(il)
```

The immunity values should represent the proportion of the population
with sufficient neutralising antibody titres to prevent infection by
each lineage. These are necessarily approximate — no surveillance system
measures this directly — but even rough estimates improve the
decomposition substantially over ignoring immunity entirely.

## Fitness decomposition

Given a fitted frequency model and an immunity landscape,
[`fitness_decomposition()`](https://CuiweiG.github.io/lineagefreq/reference/fitness_decomposition.md)
partitions each lineage’s growth advantage into intrinsic and escape
components:

``` r
data(cdc_sarscov2_jn1)
x <- lfq_data(cdc_sarscov2_jn1, lineage = lineage,
              date = date, count = count)
fit <- fit_model(x, engine = "mlr")

fd <- fitness_decomposition(fit, il, generation_time = 5)
fd
plot(fd)
```

The decomposition follows the framework of Figgins and Bedford (2025).
For each lineage, the effective fitness is modelled as the product of
intrinsic transmissibility and the susceptible fraction of the
population. The escape component is the differential susceptibility in
log space: a variant facing lower population immunity has a fitness
advantage from escape, regardless of its intrinsic transmissibility.

## DMS-informed early detection

When a new variant is first detected, sequence counts are too sparse for
reliable growth rate estimation. Deep Mutational Scanning (DMS) data —
laboratory measurements of how individual mutations affect receptor
binding and antibody escape — can serve as an informative prior.
[`fit_dms_prior()`](https://CuiweiG.github.io/lineagefreq/reference/fit_dms_prior.md)
implements penalised multinomial logistic regression that shrinks growth
rate estimates toward DMS-derived priors:

``` r
# DMS escape scores from Bloom Lab measurements
dms <- c("JN.1" = 0.05, "XBB.1.5" = 0.01, "BA.5" = -0.02)

fit_dms <- fit_dms_prior(x, dms_scores = dms, lambda = 3)
growth_advantage(fit_dms)
```

The `lambda` parameter controls the strength of the prior. At
`lambda = 0`, the result is identical to standard MLR. As `lambda`
increases, estimates are pulled toward the DMS prior. The appropriate
value depends on the relative informativeness of the DMS data and the
quantity of observed sequences — a judgment call that should be informed
by cross-validation or the calibration tools from the previous vignette.

## Selective pressure as an early warning signal

The
[`selective_pressure()`](https://CuiweiG.github.io/lineagefreq/reference/selective_pressure.md)
function computes a population-level metric from genomic data alone,
requiring no case counts or epidemiological data:

``` r
sp <- selective_pressure(fit)
sp

# The variance method captures selection intensity
sp_var <- selective_pressure(fit, method = "variance")
```

When selective pressure is positive and increasing, the
population-average fitness of circulating viruses is rising. This
precedes epidemic growth because it means the effective reproduction
number of the average virus is increasing. The metric is particularly
valuable in settings where case reporting is incomplete or delayed — a
common situation in low-resource surveillance systems.

## Practical considerations

The fitness decomposition is only as reliable as the immunity estimates
that feed into it. Uncertainty in seroprevalence data propagates
directly into uncertainty about the transmissibility versus escape
partition. We recommend treating the decomposition as a structured
hypothesis about the mechanisms driving variant replacement, not as a
precise measurement.

DMS priors work best when the mutations in a new variant have been
characterised in the laboratory. For variants with entirely novel
mutations, the DMS prior carries no information and the method reduces
to standard MLR.

The selective pressure metric assumes that the MLR model adequately
captures the variant fitness landscape. During periods of rapid immune
landscape change (e.g., mass vaccination campaigns), the
constant-growth-rate assumption may be violated, and the metric should
be interpreted cautiously.

## References

- Figgins MD, Bedford T (2025). Jointly modeling variant frequencies and
  case counts. *medRxiv*. <doi:10.1101/2024.12.02.24318334>
- Dadonaite B et al. (2023). Deep mutational scanning of the full
  SARS-CoV-2 spike. *Cell* 186(6):1263–1278.
- Bloom JD, Neher RA (2023). Fitness effects of mutations to SARS-CoV-2
  proteins. *Virus Evolution* 9(2):vead055.
