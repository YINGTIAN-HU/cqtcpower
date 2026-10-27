# C-QTc Power and Sample Size for Parallel SAD Studies

R code for evaluating **power** and **sample size** of concentration–QTc (C-QTc)
analyses in parallel-group single ascending dose (SAD) studies, using linear
mixed-effects models (LMMs).

The project does two things:

1. **Analytical power / sample size.** Closed-form expressions for the variance
   of the fixed effects (theta_Z, theta_C) in random-intercept and random-slope
   LMMs, turned into power and sample-size formulas.
2. **Simulation validation.** Simulate PK concentrations and ΔQTc data, fit the
   LMMs with `lme4`/`lmerTest`, and compare empirical power, coverage
   probability (CP), and relative bias with the analytical results.

The quantity being tested is the predicted ΔΔQTc at the supratherapeutic
concentration `C_supra`:

    U = theta_Z + theta_C * C_supra       (models with a treatment effect Z)
    U = theta_C * C_supra                 (models without Z)

A study "succeeds" (no QT effect) when the upper bound of the one-sided 95% CI
is below `U0 = 10` ms.

---

## Folder layout

```
├── README                         this file
│
├── utils.R                        low-level matrix helpers (gi, ai, bi, A^-1, ...)
├── generate_concentration.R       PK concentration generators
├── generate_delta_QTC.R           ΔQTc generator for one cohort
├── generate_parallel_SAD_data.R   full parallel SAD study + LMM fit
├── calculate_variance.R           analytical (inverse) covariance of (theta_Z, theta_C)
├── calculate_power.R              analytical + empirical power, CP, bias
├── calculate_analytical_sample_size.R   analytical sample size
│
├── simulation_power_evaluation_cluster_script.R   one simulation scenario (run on cluster)
├── aggregate_simulation_results.R  combines all scenarios into Results/power.txt
├── table_results_summary.R         prints Tables 4, 5, 6
├── plot_simulate_QTc.R             example simulated data figure (Figure/QTC.pdf)
├── plot_sample_size_estimation_comparison.R   sample-size figures (Figures 2–4)
├── plot_figure_appendix.R          appendix power figures (Figures B1, B2)
│
├── Results/                       simulation output, one folder per scenario
├── Figure/                        generated PDF figures
```

---

## How the code fits together

```
                 utils.R
                    │
generate_concentration.R ──► generate_delta_QTC.R ──► generate_parallel_SAD_data.R
                                                         │ (simulate + fit LMM)
                    calculate_variance.R ◄── utils.R     │
                         │                               │
          ┌──────────────┴──────────────┐                │
calculate_analytical_sample_size.R   calculate_power.R ◄─┘
          │                               │
plot_sample_size_estimation_        simulation_power_evaluation_cluster_script.R
  comparison.R  (Figs 2–4)                │  writes Results/<scenario>/{power,CP,bias}.txt
                                          ▼
                       aggregate_simulation_results.R ──► Results/power.txt
                                          │
                  table_results_summary.R (Tables 4–6), plot_figure_appendix.R (Figs B1–B2)
```

---

## Models

Six LMMs are compared throughout (names used as `model =` arguments and as row
labels in the results files). Every model also has fixed effects for baseline
QTc (`Baseline`) and time (`t1 … tT`), and no global intercept.

| Label     | Fixed effects            | Random effects (by subject)              |
|-----------|--------------------------|------------------------------------------|
| `RI`      | Concentration            | intercept                                |
| `RI + Z`  | Treatment + Concentration| intercept                                |
| `RS0`     | Concentration            | intercept + slope, uncorrelated (`||`)   |
| `RS0 + Z` | Treatment + Concentration| intercept + slope, uncorrelated (`||`)   |
| `RS`      | Concentration            | intercept + slope, correlated            |
| `RS + Z`  | Treatment + Concentration| intercept + slope, correlated            |

Analytical power is available for the first four only; `RS` / `RS + Z` are
`NA` in the analytical column.

---

## Getting started

### Requirements
R (≥ 4.0) with the following packages:

```r
install.packages(c("lme4", "lmerTest", "nlme", "MASS",
                   "dplyr", "tidyr", "ggplot2", "cowplot", "scales",
                   "gridExtra", "latex2exp"))
```

`pbkrtest` is also needed for Kenward–Roger df in `lmerTest::contest`.

### Typical workflows
Open `Bala.Rproj` first so the working directory is the project root.

**1. Analytical sample size (fast, no simulation)**
```r
source("utils.R"); source("calculate_variance.R"); source("calculate_analytical_sample_size.R")
calculate_analytical_sample_size2(N_G = 8, R = 3, delta = 5, sigma_e = 6, sigma0 = 3, sigma1 = 600)
```

**2. Power for one scenario (analytical + simulation)**
```r
library(lme4); library(lmerTest); library(MASS)
for (f in c("utils.R", "generate_concentration.R", "generate_delta_QTC.R",
            "generate_parallel_SAD_data.R", "calculate_variance.R", "calculate_power.R"))
  source(f)
res <- calculate_power(U1 = 5, CV = 0.75, theta_Z_true = 0.5, N_sim = 100)
res$power; res$CP; res$bias
```
Each empirical fit is an `lmer` call, and the code runs `6 × N_sim` fits. Start
with a small `N_sim`.

**3. Reproduce the manuscript results**
1. Run `simulation_power_evaluation_cluster_script.R` for every scenario on the cluster.
2. Copy the scenario folders into `Results/`.
3. Run `aggregate_simulation_results.R` → `Results/power.txt`.
4. Run `table_results_summary.R` (Tables 4–6) and `plot_figure_appendix.R` (Figures B1–B2).
5. Run `plot_sample_size_estimation_comparison.R` (Figures 2–4) and `plot_simulate_QTc.R`.

---

## Default simulation parameters

| Parameter | Default | Meaning |
|---|---|---|
| `dose_level` | 1, 3, 10, 30, 100, 200, 400, 800 mg | 8 SAD cohorts |
| `N_A`, `N_P` | 6, 2 | active / placebo subjects per cohort |
| `t` | 1, 2, 4, 8, 12, 16, 24 h | ECG time points |
| `Ka, CL, V` | 0.8, 12000, 40000 | one-compartment PK parameters |
| `theta0, theta_W` | -5.76, -0.58 | intercept, baseline effect |
| `theta_T` | 1.93, 2.29, 3.53, -3.15, -0.28, 5.17, 3.98 | time effects |
| `eta0` | 8 | SD of random intercept |
| `epsilon` | 5 | residual SD |
| `mean_baseline`, `sd_baseline` | 408, 16 | baseline QTc |
| `C_supra` | 0.008 (also 0.0113) | supratherapeutic concentration |
| `U0` | 10 ms | regulatory threshold |

---
