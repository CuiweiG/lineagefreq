# Twitter/X Announcement Thread — lineagefreq

**Tweet 1:** lineagefreq is now on CRAN. An R package for estimating
pathogen growth advantages, fitting variant frequency trajectories, and
forecasting lineage replacement dynamics from genomic surveillance
counts. <https://CRAN.R-project.org/package=lineagefreq> \#rstats
\#genomics \#PublicHealth

**Tweet 2:** Three lines from raw counts to growth advantage estimates:

``` r
x <- lfq_data(cdc_data, lineage, date, count)
fit <- fit_model(x, engine = "mlr")
growth_advantage(fit, type = "relative_Rt", generation_time = 5)
```

Ships with real CDC data for immediate reproducibility.

**Tweet 3:** Five engines behind one interface — frequentist MLR,
hierarchical pooling for sparse multi-site data, Piantham Rt conversion,
and two Bayesian models via Stan. Rolling-origin backtesting built in so
you evaluate forecasts honestly, not on training data. \#OpenSource
\#Epidemiology

**Tweet 4:** Validated against published estimates: BA.2 vs BA.1
relative Rt ≈ 1.34 (cf. Lyngse et al. 2022). 2-week forecast MAE ~4
pp. Includes sequencing power calculations and emerging lineage
detection. Four vignettes with real-world case studies.

**Tweet 5:** install.packages(“lineagefreq”) Built for surveillance
teams, epi labs, and public health analysts who need reproducible
variant tracking without assembling bespoke pipelines. — @CuiweiG23
\#GenomicSurveillance \#SARSCoV2 \#rstats
