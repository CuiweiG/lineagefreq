# Joint Conformal Prediction on the Simplex

Produces a compositional prediction region that respects the simplex
constraint (frequencies sum to 1) using Aitchison geometry. Unlike
marginal conformal prediction, joint prediction guarantees that the
*entire* frequency vector is covered, not just individual lineages.

## Usage

``` r
conformal_forecast_joint(
  fit,
  data,
  horizon = 28L,
  ci_level = 0.95,
  cal_fraction = 0.3,
  seed = NULL
)
```

## Arguments

- fit:

  An `lfq_fit` object.

- data:

  An `lfq_data` object (same data used to fit the model).

- horizon:

  Integer; forecast horizon in days (default 28).

- ci_level:

  Numeric in (0,1); coverage target (default 0.95).

- cal_fraction:

  Numeric in (0,1); fraction of dates for calibration (default 0.3).

- seed:

  Optional integer for reproducibility.

## Value

An `lfq_conformal_joint` S3 object (list) with:

- forecast:

  The point forecast (`lfq_forecast`).

- radius:

  Conformal radius in Aitchison distance.

- marginal_intervals:

  Tibble with .date, .lineage, .lower_joint, .upper_joint — marginal
  bounds projected from the joint region.

- marginal_only:

  Tibble with .lower_marginal, .upper_marginal from standard marginal
  conformal prediction.

- comparison:

  Tibble indicating whether joint intervals are wider or narrower than
  marginal intervals per lineage.

- calibration_scores:

  Aitchison distances on calibration set.

- ci_level:

  Nominal coverage level.

- n_cal:

  Number of calibration compositions.

## Details

The nonconformity score is the Aitchison distance between predicted and
observed compositions. The Aitchison distance equals the Euclidean
distance in the isometric log-ratio (ILR) transformed space, which
respects the geometry of the simplex (Aitchison, 1986).

The prediction region is the set of all compositions within Aitchison
distance \\r\\ of the point forecast, where \\r\\ is the
\\(1-\alpha)(1+1/n)\\ quantile of calibration distances. Marginal
intervals are obtained by projecting this region onto each coordinate
axis, then intersecting with \\\[0, 1\]\\.

## References

Aitchison J (1986). *The Statistical Analysis of Compositional Data*.
Chapman & Hall.

Vovk V, Gammerman A, Shafer G (2005). *Algorithmic Learning in a Random
World*. Springer.
