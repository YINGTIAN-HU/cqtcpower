# cqtcpower

Power and sample size for concentration-QTc (C-QTc) analyses of parallel-group
single ascending dose (SAD) studies analysed with linear mixed-effects models.

A study concludes **no clinically relevant QT effect** when the upper bound of
the one-sided 95% confidence interval of the predicted ΔΔQTc at the
supratherapeutic concentration `C_supra` is below 10 ms. `cqtcpower` answers:

* **Power:** how likely is that conclusion for a given number of subjects per
  cohort? Use `cqtc_power()` and `cqtc_power_curve()`.
* **Sample size:** how many subjects per cohort reach a target power? Use
  `cqtc_sample_size()`.

## Installation

```r
# from the folder that contains cqtcpower/
install.packages("cqtcpower", repos = NULL, type = "source")
```

Requires `lme4`, `lmerTest` and `MASS`. Kenward-Roger degrees of freedom need
`pbkrtest`; parallel simulation needs `future.apply`.

## Quick start

```r
library(cqtcpower)

des <- cqtc_design(
  conc = conc_pk_1cmt(dose = c(1, 3, 10, 30, 100, 200, 400, 800)),
  n_active = 6, n_placebo = 2,
  times = c(1, 2, 4, 8, 12, 16, 24)
)
par <- cqtc_params(U1 = 5, C_supra = 0.008, theta_Z = 0.5,
                   sigma_e = 5, sigma0 = 8, slope_cv = 0.75)

cqtc_power(des, par)                  # power with 6 + 2 subjects per cohort
cqtc_sample_size(des, par, power = 0.8)
plot(cqtc_power_curve(des, par))

# simulation check, including the correlated random-slope models
sim <- cqtc_simulate(des, par, nsim = 1000, ddf = "Kenward-Roger")
```

Output of `cqtc_sample_size(des, par, power = 0.8)`:

```
C-QTc sample size (analytical), target power 80%
Design: 8 cohorts, active:placebo = 6:2, 7 ECG times
Truth:  U1 = 5 ms at C_supra = 0.008 (theta_Z = 0.5); sigma_e = 5, sigma0 = 8, sigma1 = 421.9
Rule:   one-sided 95% upper bound < U0 = 10 ms

   model n/cohort (A:P) total N power SE(U)
      RI        4 (3:1)      32 0.819 1.956
  RI + Z      16 (12:4)     128 0.847 1.873
     RS0        8 (6:2)      64 0.888 1.748
 RS0 + Z      20 (15:5)     160 0.870 1.803
```

## Interactive app

```r
install.packages(c("shiny", "bslib"))   # once
cqtcpower::cqtc_app()
```

Set the design and assumptions in the sidebar. The tabs show power at the
current cohort size and a power curve, the sample size for the target power
(also across true ΔΔQTc values), a simulation check of any of the six models,
and the typical concentration profiles against `C_supra`. Tables can be
downloaded as CSV.

## Main functions

| Function | Purpose |
|---|---|
| `conc_pk_1cmt()`, `conc_profile()` | Concentration model (one-compartment PK, or a fixed profile scaled per cohort) |
| `cqtc_design()` | Cohorts, active/placebo allocation, ECG times |
| `cqtc_params()` | True ΔΔQTc `U1`, treatment effect, variance components, threshold and α |
| `cqtc_power()` | Analytical power for `RI`, `RI + Z`, `RS0`, `RS0 + Z` |
| `cqtc_power_curve()` | Power across cohort sizes, with `plot()` |
| `cqtc_sample_size()` | Smallest cohort size reaching the target power (analytical or simulation) |
| `cqtc_simulate()` | Empirical power, coverage and bias for all six models (Wald / Satterthwaite / Kenward-Roger) |
| `simulate_sad_data()`, `fit_cqtc_model()` | Simulate one study; fit one model and get the bound at `C_supra` |

## Models

| Model | Fixed effects | Random effects (by subject) | Analytical |
|---|---|---|---|
| `RI` | concentration | intercept | yes |
| `RI + Z` | treatment + concentration | intercept | yes |
| `RS0` | concentration | intercept + slope, uncorrelated | yes |
| `RS0 + Z` | treatment + concentration | intercept + slope, uncorrelated | yes |
| `RS` | concentration | intercept + slope, correlated | simulation only |
| `RS + Z` | treatment + concentration | intercept + slope, correlated | simulation only |

All models also include fixed effects for baseline QTc and each time point.

## Differences from the original research code

* Random-slope sample size uses C* = C / (1 + g_i), as derived in Appendix A.2
  (the original `calculate_analytical_sample_size2()` used `sqrt(1 + g_i)`).
* The term g_i²/(1 + g_i)² is averaged over all simulated concentration
  datasets (the original kept only the last one).
* Analytical and simulation results use one set of model names.
* All models are fitted to the same simulated dataset in each replicate, which
  is equivalent to the original per-model seeding but six times cheaper.
* Models are fitted on `C_rel = Concentration / C_supra` (an exact
  reparameterisation), with fallback optimisers. In the paper's design this
  removes the convergence warnings that ~8% of `RS0` and ~18% of `RS` fits
  produced on the mg/mL scale. Convergence and boundary (singular) fits are
  reported per model by `cqtc_simulate()`.
