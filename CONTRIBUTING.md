# Contributing to lineagefreq

Thank you for your interest in contributing to lineagefreq! This
document outlines how to propose changes and submit contributions.

## Reporting Bugs

If you find a bug, please [open an
issue](https://github.com/CuiweiG/lineagefreq/issues/new) with a minimal
reproducible example. The [reprex](https://reprex.tidyverse.org/)
package is the preferred way to create one. Include:

- A brief description of the expected vs. actual behaviour
- A reprex that demonstrates the problem
- Output of [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
  or `devtools::session_info()`

## Suggesting Features

Feature requests are welcome. Please open an issue describing the use
case and any alternatives you have considered.

## Pull Requests

To contribute code:

1.  Fork the repository and create a new branch from `main`.
2.  If adding a new feature, include tests in `tests/testthat/`.
3.  Make sure `devtools::check()` passes with no errors or warnings.
4.  Use the [tidyverse style guide](https://style.tidyverse.org/). You
    can check your code with `styler::style_pkg()`.
5.  Update documentation with `devtools::document()` if you change any
    roxygen comments.
6.  Submit a pull request describing the changes.

## Code Style

This package follows the [tidyverse style
guide](https://style.tidyverse.org/). Please use
[styler](https://styler.r-lib.org/) to format code and
[lintr](https://lintr.r-lib.org/) to check for common issues.

## Testing

All new functionality should be accompanied by
[testthat](https://testthat.r-lib.org/) tests. Run the test suite with:

``` r
devtools::test()
```

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://cuiweig.github.io/lineagefreq/CODE_OF_CONDUCT.md). By
participating you agree to abide by its terms.
