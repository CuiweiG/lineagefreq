# Validation Analysis

Supplementary analysis scripts for the lineagefreq R package.
These scripts reproduce all validation results reported in the
manuscript.

## Requirements

- R >= 4.1.0
- lineagefreq (development version from this repository)
- Additional packages: future, furrr, viridis, kableExtra,
  gridExtra, scales

Install with:
```r
install.packages(c("future", "furrr", "viridis", "kableExtra",
                   "gridExtra", "scales"))
```

## Running

From the package root directory:

```r
setwd("path/to/lineagefreq")
source("analysis/run_validation.R")
```

Or from the command line:

```bash
cd lineagefreq
Rscript analysis/run_validation.R
```

## Expected runtime

| System | Cores | Estimated time |
|--------|-------|---------------|
| AMD EPYC 9654 (192 cores) | 172 | ~5 minutes |
| Standard workstation (8 cores) | 7 | ~15 minutes |
| Laptop (4 cores) | 3 | ~25 minutes |

## Output

```
analysis/
├── results/
│   ├── benchmark_raw.rds      # Backtest results all engines
│   ├── calibration_raw.rds    # PIT values and reliability data
│   ├── surveillance_raw.rds   # EVOI and allocation results
│   └── decomposition_raw.rds  # Fitness decomposition
├── figures/
│   ├── fig1_benchmark_heatmap.pdf
│   ├── fig2_benchmark_errorbar.pdf
│   ├── fig3_pit_histogram.pdf
│   ├── fig4_reliability_diagram.pdf
│   ├── fig5_conformal_vs_parametric.pdf
│   ├── fig6_adaptive_allocation.pdf
│   ├── fig7_alert_threshold.pdf
│   ├── fig8_fitness_decomposition.pdf
│   ├── fig9_influenza_demo.pdf
│   └── fig10_graphical_abstract.pdf
└── tables/
    ├── table1_benchmark.tex
    ├── table2_calibration.tex
    └── table3_surveillance.tex
```

## Data sources

All analyses use datasets shipped with the lineagefreq package:

- `cdc_ba2_transition`: Real CDC SARS-CoV-2 variant proportions,
  BA.1 to BA.2 transition (Dec 2021 – Jun 2022). Public domain.
- `cdc_sarscov2_jn1`: Real CDC biweekly data, JN.1 emergence
  (Oct 2023 – Jun 2024). Public domain.
- `influenza_h3n2`: Simulated H3N2 clade data based on
  Nextstrain dynamics.
