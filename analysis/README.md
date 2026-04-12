# lineagefreq Validation Analysis — Nature Methods

## Execution

Run the full pipeline from the package root directory:

```r
source("analysis/run_all.R")
```

Or run individual scripts in order:

```r
source("analysis/00_setup.R")
source("analysis/01_data_prep.R")
# ... etc.
```

**Important:** Scripts must be run in numerical order. Each depends on outputs from prior scripts.

## Expected Runtime (AMD EPYC 9654, 96 cores, 192 GB RAM)

| Script | Description | Est. Time |
|--------|------------|-----------|
| 00_setup.R | System config, theme | <5 sec |
| 01_data_prep.R | All data preprocessing | 1-2 min |
| 02_benchmark.R | Multi-country backtest | 5-15 min |
| 03_calibration.R | PIT, reliability, conformal | 3-5 min |
| 04_decision_impact.R | Trigger analysis, sample size, window | 5-10 min |
| 05_fitness.R | Fitness decomposition + sensitivity | 2-3 min |
| 06_surveillance.R | Adaptive allocation simulation | 3-5 min |
| 07_influenza.R | Multi-pathogen demonstration | 1-2 min |
| 08_figures.R | Nature Methods 3-figure layout | 1-2 min |
| 09_tables.R | LaTeX booktabs tables | <30 sec |
| **Total** | | **~20-40 min** |

## Output Manifest

### Results (analysis/results/)

| File | Description |
|------|-------------|
| ecdc_prepared.rds | ECDC multi-country lfq_data objects |
| cdc_regional_prepared.rds | CDC HHS regional lfq_data objects |
| immune_landscape.rds | Time-varying immunity estimates |
| flunet_prepared.rds | WHO FluNet influenza data |
| benchmark_multicountry.rds | Backtest metrics, runtime, Bedford comparison |
| calibration_comparison.rds | PIT, reliability, conformal/recalibrated |
| decision_impact.rds | Vaccine trigger timing analysis |
| sample_size_analysis.rds | MAE and coverage vs sample size |
| window_analysis.rds | Training window vs PIT uniformity |
| fitness_sensitivity.rds | Fitness decomposition + perturbation |
| surveillance_simulation.rds | Adaptive allocation simulation results |
| evoi_results.rds | Expected Value of Information curve |
| influenza_results.rds | Multi-pathogen demonstration results |
| alert_threshold.rds | SPRT alert retrospective |

### Figures (analysis/figures/)

| File | Description | Dimensions |
|------|-------------|-----------|
| figure1.pdf / .png | Overview + Benchmark (6 panels) | 180 × 220 mm |
| figure2.pdf / .png | Calibration (6 panels) | 180 × 250 mm |
| figure3.pdf / .png | Advanced analyses (5 panels) | 180 × 200 mm |

### Tables (analysis/tables/)

| File | Description |
|------|-------------|
| table1.tex | Multi-country benchmark with Bedford Lab reference |
| table2.tex | Calibration comparison (Parametric/Conformal/Recalibrated) |
| table3.tex | Decision impact: vaccine trigger timing |
| table4.tex | Sample size requirements for MAE and coverage |

## Data Provenance

| Dataset | Source | URL | License |
|---------|--------|-----|---------|
| cdc_variant_proportions_full.csv | CDC COVID Data Tracker, Variant Proportions | https://data.cdc.gov/Laboratory-Surveillance/SARS-CoV-2-Variant-Proportions/jr58-6ysp | Public Domain (US Government) |
| ecdc_variants.csv | ECDC TESSy, Genomic surveillance | https://www.ecdc.europa.eu/en/publications-data/data-virus-variants-covid-19-eueea | ECDC Copyright, reuse permitted |
| who_flunet.csv | WHO FluNet | https://www.who.int/tools/flunet | WHO Terms of Use |
| owid_vaccinations.csv | Our World in Data (OWID) | https://github.com/owid/covid-19-data | CC BY 4.0 |
| cdc_seroprevalence.csv | CDC Nationwide Blood Donor Seroprevalence Survey | https://covid.cdc.gov/covid-data-tracker/#nationwide-blood-donor-seroprevalence | Public Domain (US Government) |
| cdc_sarscov2_ba2, cdc_sarscov2_jn1 | Package built-in | lineagefreq::cdc_sarscov2_ba2 | Per package license |

## Known Limitations

1. **Reconstructed counts (CDC regional):** Proportions converted to counts assuming 800 sequences per region per biweekly period. This is a conservative assumption; actual sequencing volume varied.

2. **FluNet data format:** WHO FluNet data may be downloaded as HTML tables rather than CSV. If parsing fails, the pipeline falls back to simulated influenza data (clearly labeled).

3. **Piantham engine:** Requires CmdStan. If unavailable, only the MLR engine is benchmarked (reported honestly as a 1-row heatmap per dataset).

4. **Immune landscape construction:** Vaccine efficacy waning model uses simplified exponential decay (half-life = 6 months). Real-world waning is heterogeneous by age, vaccine type, and prior infection history.

5. **Fitness decomposition identifiability:** The transmissibility/escape decomposition requires time-varying immunity (>5 pp change during observation window). Sensitivity analysis (±30% perturbation) quantifies this limitation.

6. **ACI coverage guarantee:** Adaptive conformal inference provides asymptotic coverage guarantees (Gibbs & Candès 2021). Finite-sample coverage may deviate under strong distribution shift.

7. **Sample size analysis:** Uses multinomial resampling from observed proportions, which assumes the original data represents "truth." This understates uncertainty from systematic biases in genomic surveillance.
