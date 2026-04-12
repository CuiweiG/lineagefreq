# Plot calibration diagnostics

Produces either a reliability diagram (default) or a PIT histogram.

## Usage

``` r
# S3 method for class 'calibration_report'
plot(x, type = c("reliability", "pit"), ...)
```

## Arguments

- x:

  A `calibration_report` object.

- type:

  Which panel to display: `"reliability"` (default) or `"pit"` for the
  PIT histogram.

- ...:

  Unused.

## Value

A ggplot object.
