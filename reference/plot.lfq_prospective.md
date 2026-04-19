# Plot Prospective Evaluation Results

Plot Prospective Evaluation Results

## Usage

``` r
# S3 method for class 'lfq_prospective'
plot(x, type = c("coverage", "radius", "comparison"), ...)
```

## Arguments

- x:

  An `lfq_prospective` object.

- type:

  Character; one of `"coverage"` (cumulative coverage over time),
  `"radius"` (conformal radius over time), or `"comparison"`
  (retrospective vs prospective).

- ...:

  Ignored.

## Value

A ggplot object.
