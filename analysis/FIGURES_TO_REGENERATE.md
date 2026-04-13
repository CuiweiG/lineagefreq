# Figures Requiring Regeneration

After the audit fixes, the following figures need to be regenerated
in RStudio. Run from the package root directory.

## Figure 3 (prospective evaluation)
- **Script:** `analysis/regen_prospective.R`
- **Fix:** Panel c x-axis labels shortened to "ACI", "Param", "Static CF"
  with 45-degree rotation
- **Run:** `source("analysis/regen_prospective.R")`
- **Output:** `submission/figures/figure_prospective.pdf`

## Figure 6 (joint conformal)
- **Script:** `analysis/make_figure_joint.R`
- **Fix:** Legend now uses distinct colours (blue/orange) for Marginal/Joint
  instead of identical black circles
- **Run:** `source("analysis/make_figure_joint.R")`
- **Output:** `submission/figures/fig06_joint_conformal.pdf`

## Notes
- All scripts require `Sys.setlocale("LC_TIME", "English")` which is
  already included
- lineagefreq >= 0.6.0 must be installed: `devtools::install(".")`
- Figure 1, engine comparison, and extended data figures do NOT need
  regeneration — they passed the audit
