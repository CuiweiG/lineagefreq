# Validation Analysis

Supplementary analysis scripts for the lineagefreq R package.
Produces all figures, tables, and benchmark results reported in
the manuscript.

## Execution order

```
00_setup.R       → Environment configuration, parallel backend
01_benchmark.R   → Rolling-origin backtest, forecast accuracy
02_calibration.R → PIT diagnostics, conformal comparison
03_surveillance.R → EVOI, adaptive allocation, alert threshold
04_fitness.R     → Immune landscape, fitness decomposition
05_influenza.R   → Multi-pathogen demonstration
06_figures.R     → Generate all 10 publication figures
07_tables.R      → Generate all 3 LaTeX tables
```

## Running

**Full analysis** (recommended):
```r
setwd("path/to/lineagefreq")
source("analysis/run_all.R")
```

**Individual sections** (after 00_setup.R):
```r
source("analysis/00_setup.R")
source("analysis/02_calibration.R")  # loads saved results from 01
```

Each script after 01 loads results from `analysis/results/`, so
you can rerun individual sections without repeating the full
backtest computation.

## Expected runtime

| System | Cores | Approx. time |
|--------|-------|-------------|
| AMD EPYC 9654 | 192 | ~2 min |
| Standard workstation | 8 | ~10 min |
| Laptop | 4 | ~20 min |

Most time is spent in 01_benchmark.R (backtesting). Sections
02-07 load saved results and run in seconds.

## Required packages

Core (from lineagefreq):
- R >= 4.1.0, devtools

Analysis-specific:
- ggplot2, dplyr, tidyr, tibble
- viridis, scales, gridExtra
- future, furrr (parallel computation)
- kableExtra (LaTeX tables)
- desc (version detection)

Install all:
```r
install.packages(c("devtools", "ggplot2", "dplyr", "tidyr",
  "tibble", "viridis", "scales", "gridExtra", "future", "furrr",
  "kableExtra", "desc"))
```

## Known limitations

- **Bayesian engines (fga, garw)** require CmdStan. If not
  installed, 01_benchmark.R benchmarks only the frequentist
  engines (mlr, hier_mlr, piantham). This is documented honestly
  in the output.
- **SPRT alerting** may not trigger on biweekly CDC data due to
  insufficient observations. This is a real limitation of
  sequential testing on low-frequency data, documented in
  03_surveillance.R.
- **Fitness decomposition** requires external immunity estimates.
  The values used in 04_fitness.R are approximate and should be
  treated as illustrative.

## Output manifest

```
analysis/
├── results/
│   ├── benchmark_ba2.rds       Backtest + scores, BA.2 dataset
│   ├── benchmark_jn1.rds       Backtest + scores, JN.1 dataset
│   ├── calibration_results.rds PIT, reliability, conformal data
│   ├── surveillance_results.rds EVOI, allocation, alerts
│   ├── fitness_results.rds     Decomposition results
│   └── influenza_results.rds   Influenza demo data
├── figures/
│   ├── fig01_benchmark.pdf     Accuracy heatmap
│   ├── fig02_mae_horizon.pdf   MAE by horizon with CI
│   ├── fig03_pit.pdf           PIT histograms (KEY FIGURE)
│   ├── fig04_reliability.pdf   Reliability diagram
│   ├── fig05_conformal.pdf     Conformal vs parametric
│   ├── fig06_surveillance.pdf  EVOI + allocation
│   ├── fig07_alert.pdf         JN.1 detection
│   ├── fig08_fitness.pdf       Fitness decomposition
│   ├── fig09_influenza.pdf     Influenza demo
│   └── fig10_abstract.pdf      Graphical abstract
│   └── (+ .png versions of each)
└── tables/
    ├── table1_benchmark.tex    Accuracy comparison
    ├── table2_calibration.tex  Calibration metrics
    └── table3_surveillance.tex Allocation comparison
```
