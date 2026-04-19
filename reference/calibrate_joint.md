# Joint Calibration Diagnostics

Assess calibration of multivariate forecasts using the energy score,
joint coverage, and multivariate rank histograms.

## Usage

``` r
calibrate_joint(bt, n_bins = 10L)
```

## Arguments

- bt:

  An `lfq_backtest` object.

- n_bins:

  Integer; number of rank histogram bins (default 10).

## Value

A `joint_calibration_report` S3 object (list) with:

- energy_scores:

  Tibble with origin_date, horizon, energy_score.

- mean_energy_score:

  Overall mean energy score.

- joint_coverage:

  Tibble with nominal level, observed coverage, using Aitchison distance
  from backtest results.

- rank_histogram:

  Tibble with bin, count, density, expected.

- n:

  Number of forecast-observation vectors.

## Details

The energy score is a multivariate proper scoring rule that generalises
the CRPS (Gneiting & Raftery, 2007, Section 5). For a deterministic
forecast \\f\\ and observation \\y\\: \$\$ES = \|\|f - y\|\|^2\$\$ where
the norm is the Aitchison distance on the simplex.

Joint coverage at level \\1-\alpha\\: the fraction of observed
composition vectors falling within Aitchison distance \\q\_\alpha\\ of
the forecast, where \\q\_\alpha\\ is derived from the empirical
distribution of distances.

The multivariate rank histogram uses the pre-rank approach: the rank of
the observation among an ensemble is computed using the Aitchison
distance to the ensemble mean.
