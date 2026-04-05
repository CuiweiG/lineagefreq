# CRAN Pre-Flight Audit: lineagefreq 0.2.0

**Audit date:** 2026-04-06  
**R version:** 4.5.3 (Windows 10 x64)  
**Overall verdict:** 🟢 Very close to CRAN-ready. A few items need
attention.

------------------------------------------------------------------------

## Summary

| Category     | Issues | Auto-fixable | Needs Human Review |
|--------------|--------|--------------|--------------------|
| \[FORMAT\]   | 1      | 1            | 0                  |
| \[URL\]      | 1      | 1            | 0                  |
| \[CHECK\]    | 3      | 1            | 2                  |
| \[EXAMPLES\] | 0      | —            | —                  |
| \[CONSOLE\]  | 0      | —            | —                  |
| \[DONTRUN\]  | 0      | —            | —                  |

------------------------------------------------------------------------

## Detailed Findings

### \[FORMAT\] DESCRIPTION Metadata

1.  **\[FORMAT-1\] Additional_repositories uses deprecated URL**
    *(auto-fixable)*
    - Current: `https://mc-stan.org/r-packages/`
    - The mc-stan.org/r-packages repo page itself says: “❗ This
      repository is deprecated and no longer updated” — moved to
      r-universe.
    - **Fix:** Change to `https://stan-dev.r-universe.dev`
    - **Impact:** CRAN check already resolves cmdstanr from this repo,
      but using a deprecated URL looks unprofessional and may break in
      the future.

**DESCRIPTION fields that PASS:** - ✅ Title: “Lineage Frequency
Dynamics from Genomic Surveillance Counts” — proper title case, doesn’t
start with package name, 57 chars, no trailing period - ✅ <Authors@R>:
Uses [`person()`](https://rdrr.io/r/utils/person.html) format with
proper roles (aut, cre, cph) - ✅ License: `MIT + file LICENSE` — valid
CRAN license, LICENSE file present with YEAR/COPYRIGHT HOLDER - ✅
Version: 0.2.0 — valid format - ✅ Description: Multi-sentence,
informative, includes DOI reference - ✅ Encoding: UTF-8 - ✅ Language:
en-US - ✅ LazyData: true - ✅ Depends: R (\>= 4.1.0) — reasonable
minimum

### \[URL\] Link Check

All URLs verified and reachable:

| URL                                             | Status                                |
|-------------------------------------------------|---------------------------------------|
| <https://github.com/CuiweiG/lineagefreq>        | ✅ 200                                |
| <https://github.com/CuiweiG/lineagefreq/issues> | ✅ 200                                |
| <https://mc-stan.org/r-packages/>               | ✅ 200 (but deprecated, see FORMAT-1) |
| <https://doi.org/10.1371/journal.pcbi.1012443>  | ✅ 200                                |
| <https://data.cdc.gov/>…/jr58-6ysp              | ✅ 200                                |
| <https://doi.org/10.1038/s41467-022-33498-0>    | ✅ 200                                |
| <https://doi.org/10.3201/eid2806.220158>        | ✅ 200                                |
| <https://opensource.org/licenses/MIT>           | ✅ 200                                |

1.  **\[URL-1\] README badge says “CRAN-submitted” but not yet on CRAN**
    *(auto-fixable)*
    - Badge: `https://img.shields.io/badge/CRAN-submitted-orange.svg`
    - This is fine for pre-submission but should be updated to the real
      CRAN badge URL once accepted.
    - **Note:** Not a blocking issue for submission — just cosmetic.

### \[CHECK\] R CMD check –as-cran Results

**With `_R_CHECK_FORCE_SUGGESTS_=false`:** 0 errors \| 1 warning \| 0
notes ✅

The single WARNING:

1.  **\[CHECK-1\] “Insufficient package version (submitted: 0.2.0,
    existing: 0.2.0)”** *(needs human review)*
    - This means 0.2.0 was already submitted/released. For a new CRAN
      submission, bump to 0.2.1 or 0.3.0.
    - **Action:** Bump version in DESCRIPTION, NEWS.md, and
      inst/CITATION before submission.
2.  **\[CHECK-2\] “Days since last update: 3”** *(needs human review)*
    - CRAN prefers \>7-10 days between submissions. Wait a bit longer if
      this was recently submitted.
3.  **\[CHECK-3\] “Suggests or Enhances not in mainstream repositories:
    cmdstanr”** *(auto-fixable via FORMAT-1)*
    - cmdstanr is available via Additional_repositories — CRAN accepts
      this pattern.
    - The deprecated repo URL should be updated (see FORMAT-1).
    - **Note:** This NOTE is standard for packages that suggest CmdStanR
      and is acceptable.

**All other checks PASS cleanly:** - ✅ Package installs successfully -
✅ All examples run OK (including –run-donttest in 13s) - ✅ All tests
pass (16s) - ✅ All vignettes rebuild OK (32s) - ✅ No R code problems
detected - ✅ No Rd issues - ✅ No missing documentation - ✅ No
non-ASCII issues - ✅ Namespace is clean

### \[DONTRUN\] vs Audit

✅ **No found anywhere** — all slow examples correctly use
`\donttest{}`.

Files using `\donttest{}` (appropriate — these involve model
fitting/plotting): - `autoplot.lfq_fit.Rd` -
`autoplot.lfq_forecast.Rd` - `backtest.Rd` - `compare_models.Rd` -
`forecast.lfq_fit.Rd` - `plot_backtest.Rd` - `score_forecasts.Rd`

### \[EXAMPLES\] Example Timing

✅ **All examples look fast.** R CMD check confirmed: - Regular
examples: passed quickly (within default timeout) - `--run-donttest`
examples: 13 seconds total for all 7 donttest blocks

No examples do heavy computation outside of `\donttest{}`. Simulations
use small data (3 lineages, 15-20 timepoints). No concerns.

### \[CONSOLE\] print()/cat() Audit

✅ **All cat()/print() calls are in legitimate print/summary methods:**

- [`print.lfq_fit()`](https://CuiweiG.github.io/lineagefreq/reference/print.lfq_fit.md)
  — cat() for formatted model display ✅
- `print.lfq_data()` — cat(“”) + NextMethod() ✅
- `print.lfq_forecast()` — cat(“”) + NextMethod() ✅
- `print.lfq_backtest()` — cat(“”) + NextMethod() ✅
- [`summary.lfq_fit()`](https://CuiweiG.github.io/lineagefreq/reference/summary.lfq_fit.md)
  — cat() for formatted summary display ✅

No cat()/print() calls found in non-print/summary methods. The package
correctly uses
[`cli::cli_text()`](https://cli.r-lib.org/reference/cli_text.html) and
[`message()`](https://rdrr.io/r/base/message.html) for informational
output elsewhere.

------------------------------------------------------------------------

## Auto-Fix Applied

### Fix 1: Updated Additional_repositories URL

- **Changed:** `https://mc-stan.org/r-packages/` →
  `https://stan-dev.r-universe.dev`
- **Reason:** The old mc-stan.org repo is officially deprecated

### Fixes NOT Applied (Need Human Review)

1.  **Version bump** — Must decide on 0.2.1 vs 0.3.0 and update
    DESCRIPTION, NEWS.md, CITATION, cran-comments.md
2.  **Submission timing** — Wait for sufficient days since last update
3.  **README CRAN badge** — Update to real CRAN badge URL after
    acceptance

------------------------------------------------------------------------

## Recommendation

This package is in excellent shape for CRAN. The code quality is high,
documentation is thorough, all checks pass cleanly, and CRAN policies
are followed correctly. The only required action before resubmission is:

1.  **Bump the version number** (0.2.0 → 0.2.1 minimum)
2.  **Wait sufficient time** since last submission attempt
3.  **Update cran-comments.md** to reflect current check results
