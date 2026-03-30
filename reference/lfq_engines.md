# List available modeling engines

Returns information about all modeling engines available in lineagefreq.
Core engines (mlr, hier_mlr, piantham) are always available. Bayesian
engines (fga, garw) require 'CmdStan'.

## Usage

``` r
lfq_engines()
```

## Value

A tibble with columns: `engine`, `type`, `time_varying`, `available`,
`description`.

## Examples

``` r
lfq_engines()
#> # A tibble: 5 × 5
#>   engine   type        time_varying available description                       
#>   <chr>    <chr>       <lgl>        <lgl>     <chr>                             
#> 1 mlr      frequentist FALSE        TRUE      Multinomial logistic regression (…
#> 2 hier_mlr frequentist FALSE        TRUE      Hierarchical MLR with empirical B…
#> 3 piantham frequentist FALSE        TRUE      Piantham Rt approximation from ML…
#> 4 fga      bayesian    FALSE        FALSE     Fixed growth advantage (Stan)     
#> 5 garw     bayesian    TRUE         FALSE     Growth advantage random walk (Sta…
```
