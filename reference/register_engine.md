# Register a custom modeling engine

Allows third-party packages to register custom engines with
[`fit_model()`](https://cuiweig.github.io/lineagefreq/reference/fit_model.md).
This enables an extensible plugin architecture similar to 'parsnip'
engine registration.

## Usage

``` r
register_engine(
  name,
  fit_fn,
  description = "",
  type = "frequentist",
  time_varying = FALSE
)
```

## Arguments

- name:

  Engine name (character scalar).

- fit_fn:

  Function with signature `function(data, pivot, ci_level, ...)`. Must
  return a list compatible with `lfq_fit` structure.

- description:

  One-line description of the engine.

- type:

  `"frequentist"` or `"bayesian"`.

- time_varying:

  Logical; does the engine support time-varying growth advantages?

## Value

Invisibly returns the updated engine registry.

## Examples

``` r
# Register a custom engine
my_engine <- function(data, pivot = NULL, ci_level = 0.95, ...) {
  # Custom implementation...
  .engine_mlr(data, pivot = pivot, ci_level = ci_level, ...)
}
register_engine("my_mlr", my_engine, "Custom MLR wrapper")
#> Registered engine "my_mlr".
lfq_engines()
#> # A tibble: 6 × 5
#>   engine   type        time_varying available description                       
#>   <chr>    <chr>       <lgl>        <lgl>     <chr>                             
#> 1 mlr      frequentist FALSE        TRUE      Multinomial logistic regression (…
#> 2 hier_mlr frequentist FALSE        TRUE      Hierarchical MLR with empirical B…
#> 3 piantham frequentist FALSE        TRUE      Piantham Rt approximation from ML…
#> 4 fga      bayesian    FALSE        FALSE     Fixed growth advantage (Stan)     
#> 5 garw     bayesian    TRUE         FALSE     Growth advantage random walk (Sta…
#> 6 my_mlr   frequentist FALSE        TRUE      Custom MLR wrapper                
```
