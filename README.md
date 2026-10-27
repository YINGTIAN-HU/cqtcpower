# cqtcpower

**Power and sample size for concentration-QTc (C-QTc) analyses of parallel-group
single ascending dose (SAD) studies.**

`cqtcpower` helps you plan the ECG part of a first-in-human SAD study when QT
prolongation will be assessed by C-QTc modelling with a linear mixed-effects
model. It gives you:

* **analytical power** for a given number of subjects per cohort, from closed-form
  approximations, in under a second;
* **sample size**: the smallest cohort size that reaches a target power, keeping
  your active:placebo ratio;
* **power curves** across cohort sizes;
* a **simulation engine** that fits the mixed models with `lme4`/`lmerTest` and reports
  empirical power, coverage and bias for six models (including correlated random
  slopes), with Wald, Satterthwaite or Kenward-Roger degrees of freedom;
* an interactive **Shiny app** that puts all of the above behind a point-and-click
  interface.

---

## Contents

* [Installation](#installation)
* [Background](#background)
* [Quick start](#quick-start)
* [Step by step](#step-by-step)
* [Shiny app](#shiny-app)
* [Function reference](#function-reference)
* [Methods and assumptions](#methods-and-assumptions)
* [Reproducing the paper](#reproducing-the-paper)
* [Citation](#citation) · [License](#license)

---

## Installation

```r
# from GitHub
install.packages("remotes")
remotes::install_github("<github-user>/cqtcpower")

# or from a downloaded source file
install.packages("cqtcpower_0.1.0.tar.gz", repos = NULL, type = "source")
```

Requirements: R ≥ 4.0 with `lme4`, `lmerTest` and `MASS` (installed automatically).
Optional packages:

| Package | Needed for |
|---|---|
| `pbkrtest` | Kenward-Roger degrees of freedom in simulations |
| `future.apply` | Running simulations in parallel |
| `shiny`, `bslib` | The interactive app (`cqtc_app()`) |

```r
install.packages(c("pbkrtest", "future.apply", "shiny", "bslib"))
```

---

## Background

In a parallel SAD study each cohort receives one dose level, with a few subjects
on placebo. ΔQTc is modelled against plasma concentration, and the drug is
declared to have **no clinically relevant QT effect** if the upper bound of the
one-sided 95% confidence interval of the predicted ΔΔQTc at the
supratherapeutic concentration C<sub>supra</sub> is below **U<sub>0</sub> = 10 ms**.

The quantity being tested is

* U = θ<sub>Z</sub> + θ<sub>C</sub> · C<sub>supra</sub> for models with a treatment effect Z, and
* U = θ<sub>C</sub> · C<sub>supra</sub> for models without it,

where θ<sub>C</sub> is the concentration slope and θ<sub>Z</sub> a
concentration-independent treatment effect. If the true value is U<sub>1</sub>, the
analytical power is

  power = Φ( (U<sub>0</sub> − U<sub>1</sub>) / SE(Û) − z<sub>1−α</sub> ),

with SE(Û) from closed-form approximations of the covariance of
(θ̂<sub>Z</sub>, θ̂<sub>C</sub>).

### Models

All models include fixed effects for baseline QTc and for each ECG time point
(no global intercept).

| Model | Fixed effects | Random effects (by subject) | Analytical power |
|---|---|---|---|
| `RI` | concentration | intercept | ✓ |
| `RI + Z` | treatment + concentration | intercept | ✓ |
| `RS0` | concentration | intercept + slope, uncorrelated | ✓ |
| `RS0 + Z` | treatment + concentration | intercept + slope, uncorrelated | ✓ |
| `RS` | concentration | intercept + slope, correlated | simulation only |
| `RS + Z` | treatment + concentration | intercept + slope, correlated | simulation only |

Models with Z protect against a concentration-independent drug effect at the
cost of power; random-slope models allow subjects to differ in their
concentration-QTc slope.

---

## Quick start

```r
library(cqtcpower)

# 1. the study: 8 dose cohorts, 6 active + 2 placebo each, 7 ECG times
des <- cqtc_design(
  conc = conc_pk_1cmt(dose = c(1, 3, 10, 30, 100, 200, 400, 800),
                      Ka = 0.8, CL = 12000, V = 40000, iiv_sd = 0.1),
  n_active = 6, n_placebo = 2,
  times = c(1, 2, 4, 8, 12, 16, 24)
)

# 2. what you assume to be true
par <- cqtc_params(U1 = 5, C_supra = 0.008, theta_Z = 0.5,
                   sigma_e = 5, sigma0 = 8, slope_cv = 0.75)

# 3. power with the planned design
cqtc_power(des, par)
```

```
Analytical C-QTc power
Design: 8 cohorts, active:placebo = 6:2, 7 ECG times
Truth:  U1 = 5 ms at C_supra = 0.008 (theta_Z = 0.5); sigma_e = 5, sigma0 = 8, sigma1 = 421.9
Rule:   one-sided 95% upper bound < U0 = 10 ms

   model n/cohort (A:P) total N SE(U) power
      RI        8 (6:2)      64 1.372 0.977
  RI + Z        8 (6:2)      64 2.651 0.595
     RS0        8 (6:2)      64 1.748 0.888
 RS0 + Z        8 (6:2)      64 2.853 0.543
```

```r
# 4. how many subjects per cohort for 80% power?
cqtc_sample_size(des, par, power = 0.8)
```

```
C-QTc sample size (analytical), target power 80%
...
   model n/cohort (A:P) total N power SE(U)
      RI        4 (3:1)      32 0.819 1.956
  RI + Z      16 (12:4)     128 0.847 1.873
     RS0        8 (6:2)      64 0.888 1.748
 RS0 + Z      20 (15:5)     160 0.870 1.803
```

```r
# 5. or explore everything interactively
cqtc_app()
```

---

## Step by step

### 1. Concentrations

Two ways to describe the concentrations:

```r
# one-compartment PK model with first-order absorption, one cohort per dose;
# log-normal between-subject variability on Ka, CL and V
conc_pk_1cmt(dose = c(1, 3, 10, 30, 100, 200, 400, 800),
             Ka = 0.8, CL = 12000, V = 40000, iiv_sd = 0.1)

# a fixed time profile scaled per cohort: cmax * cmax_fraction[g] * profile
conc_profile(cmax = 0.008,
             cmax_fraction = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1),
             profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01))
```

`conc_typical(design)` returns the mean profiles and the true C<sub>max</sub>,
which helps when choosing C<sub>supra</sub>. For the PK example above the top dose
reaches C<sub>max</sub> = 0.0111 at 1.96 h.

### 2. Design

```r
cqtc_design(conc, n_active = 6, n_placebo = 2, times = c(1, 2, 4, 8, 12, 16, 24))
```

When the cohort size changes, the active:placebo ratio is kept, so cohorts grow in
blocks: 4 subjects (3:1) for a 6:2 design, 3 subjects for 4:2, 6 for 5:1.

### 3. Assumptions

```r
cqtc_params(U1 = 5,            # true ddQTc at C_supra (ms)
            C_supra = 0.008,   # supratherapeutic concentration
            theta_Z = 0.5,     # concentration-independent treatment effect (ms)
            sigma_e = 5,       # residual SD (ms)
            sigma0 = 8,        # random-intercept SD (ms)
            slope_cv = 0.75,   # random-slope SD as a CV of theta_C (or give sigma1)
            U0 = 10, alpha = 0.05)
```

The slope is derived as θ<sub>C</sub> = (U<sub>1</sub> − θ<sub>Z</sub>) / C<sub>supra</sub>, and
σ<sub>1</sub> = `slope_cv` × |θ<sub>C</sub>| unless `sigma1` is given directly.

What matters for power: `U1`, `C_supra`, `theta_Z`, `sigma_e`, `sigma0`, the
random slope, `U0`, `alpha` and the design. The correlation `rho` between random
intercept and slope only enters the simulation (the analytical formulas assume
`rho = 0`). `theta0`, `theta_T`, `theta_W` and the baseline QTc distribution only
change how simulated data look; they do not affect power, so the defaults are fine.

### 4. Power and power curves

```r
cqtc_power(des, par, models = c("RI + Z", "RS0 + Z"), n_per_cohort = c(8, 16, 24))
```

```
   model n/cohort (A:P) total N SE(U) power
  RI + Z        8 (6:2)      64 2.651 0.595
  RI + Z      16 (12:4)     128 1.873 0.847
  RI + Z      24 (18:6)     192 1.529 0.948
 RS0 + Z        8 (6:2)      64 2.853 0.543
 RS0 + Z      16 (12:4)     128 2.016 0.798
 RS0 + Z      24 (18:6)     192 1.646 0.918
```

```r
pc <- cqtc_power_curve(des, par, n_per_cohort = seq(4, 40, by = 4), target = 0.8)
plot(pc)                       # power against subjects per cohort
as.data.frame(pc)              # the numbers
```

For random PK concentrations the concentration summaries are averaged over
`n_conc_sim = 1000` simulated studies (no model fitting, so this takes about a
second). Results are reproducible through `seed`.

### 5. Sample size

```r
cqtc_sample_size(des, par, power = 0.9, models = "RS0 + Z")
```

```
   model n/cohort (A:P) total N power SE(U)
 RS0 + Z      24 (18:6)     192 0.918 1.646
```

The search starts from the closed-form answer

  N = ((z<sub>1−β</sub> + z<sub>1−α</sub>) / (U<sub>0</sub> − U<sub>1</sub>))² · N<sub>ref</sub> · SE²<sub>ref</sub>(Û),

rounded up to whole blocks, and then moves block by block to the smallest cohort
size whose analytical power reaches the target. The closed-form value is kept in
the `n_formula` column.

For the correlated random-slope models, or to size on empirical power, use the
simulation search. It starts from the analytical answer and fits `nsim` studies at
each step, so it is slow:

```r
cqtc_sample_size(des, par, models = "RS + Z", method = "simulation",
                 nsim = 500, ddf = "Kenward-Roger", parallel = TRUE)
```

### 6. Simulation

```r
sim <- cqtc_simulate(des, par, nsim = 200,
                     models = c("RI + Z", "RS0 + Z", "RS + Z"),
                     ddf = c("Wald", "Satterthwaite"))
sim
```

```
Simulated C-QTc operating characteristics (200 studies, 8 subjects per cohort)
...
   model           ddf power MC SE analytical coverage rel. bias nonconv singular failed
  RI + Z          Wald 0.595 0.035      0.595    0.880     0.002       0        0      0
  RI + Z Satterthwaite 0.590 0.035      0.595    0.890     0.002       0        0      0
 RS0 + Z          Wald 0.530 0.035      0.543    0.905     0.009       0       39      0
 RS0 + Z Satterthwaite 0.515 0.035      0.543    0.910     0.009       0       39      0
  RS + Z          Wald 0.550 0.035         NA    0.900     0.003       0       48      0
  RS + Z Satterthwaite 0.540 0.035         NA    0.900     0.003       0       48      0
```

| Column | Meaning |
|---|---|
| `power` | proportion of studies whose one-sided upper bound is below U<sub>0</sub> |
| `MC SE` | Monte Carlo standard error of `power` |
| `analytical` | analytical power for comparison (where available) |
| `coverage` | proportion of two-sided (1 − 2α) intervals that contain U<sub>1</sub> |
| `rel. bias` | mean of (Û − U<sub>1</sub>) / U<sub>1</sub> |
| `nonconv` | fits with a convergence warning after trying all optimisers |
| `singular` | boundary fits (a variance component estimated as 0); valid estimates |
| `failed` | fits that stopped with an error (dropped) |

Useful options:

* `ddf`: `"Wald"` (normal quantile), `"Satterthwaite"`, `"Kenward-Roger"` (needs `pbkrtest`).
* `parallel = TRUE` with a `future` plan, e.g. `future::plan("multisession", workers = 8)`.
  Replicate *i* always uses seed `seed + i - 1`, so results do not depend on the
  number of workers.
* `n_per_cohort` to simulate a different cohort size; `drop_nonconverged = TRUE`
  to exclude non-converged fits from the summaries; `sim$replicates` holds the
  per-study estimates.

Each simulated study fits every requested model. As a guide, 1000 studies × 6
models with all three degrees-of-freedom methods took about 13 minutes on one
core (about 6.5 minutes on two) for the design above.

### 7. Lower-level functions

```r
d   <- simulate_sad_data(des, par, seed = 1)          # one simulated study (long format)
fit <- fit_cqtc_model(d, "RS0 + Z", C_supra = 0.008)   # estimate, SE and bounds at C_supra
```

`fit_cqtc_model()` can be applied to real data with the same column names
(`delta_QTC`, `Time`, `Baseline`, `Concentration`, `Treatment`, `Subject`).

---

## Shiny app

```r
cqtcpower::cqtc_app()
```

The app runs locally in your browser. Everything is set in the sidebar:

* **Study design**: one-compartment PK (doses, Ka, CL, V, between-subject SD, ECG
  times) or a fixed concentration profile; active and placebo subjects per cohort.
* **True effect and variability**: U<sub>1</sub>, C<sub>supra</sub>, θ<sub>Z</sub>, σ<sub>e</sub>,
  σ<sub>0</sub>, random slope as a CV or an SD, and ρ for simulations.
* **Decision rule**: U<sub>0</sub>, α and the target power.
* **Models** to show (RI, RI + Z, RS0, RS0 + Z).

Results update as you change the inputs. The tabs:

| Tab | What it shows |
|---|---|
| **Power** | Power of each model at the design's cohort size (green when the target is reached) and a power curve across cohort sizes. Download the curve as CSV. |
| **Sample size** | Smallest cohort size for the target power with the power achieved, and a bar chart of the required cohort size across true ΔΔQTc values. |
| **Simulation** | Runs `cqtc_simulate()` with a progress bar for any of the six models and degrees-of-freedom methods; shows empirical vs analytical power, coverage, bias and convergence counts. |
| **Concentrations** | Mean concentration-time profile of every cohort against C<sub>supra</sub>, with the true C<sub>max</sub>. |
| **About** | Decision rule, power formula and model definitions. |

**Power** — power at the planned design and the power curve:

![Power tab](man/figures/app-power.png)

**Sample size** — required cohort size for the target power, and across true ΔΔQTc:

![Sample size tab](man/figures/app-sample-size.png)

**Simulation** — empirical check by fitting the mixed models:

![Simulation tab](man/figures/app-simulation.png)

**Concentrations** — typical profiles and C<sub>supra</sub>:

![Concentrations tab](man/figures/app-concentrations.png)

Tips: power and sample size are analytical and fast; the bar chart over ΔΔQTc
takes a few seconds. Simulations fit real mixed models: 100 studies with three
models take about a minute; Kenward-Roger is slower.

---

## Function reference

| Function | Purpose |
|---|---|
| `conc_pk_1cmt()` | One-compartment PK concentration model |
| `conc_profile()` | Fixed concentration profile scaled per cohort |
| `conc_typical()` | Mean concentration profiles and true C<sub>max</sub> |
| `cqtc_design()` | Cohorts, active/placebo allocation, ECG times |
| `cqtc_params()` | True effect, variance components, threshold and α |
| `cqtc_power()` | Analytical power (RI, RI + Z, RS0, RS0 + Z) |
| `cqtc_power_curve()` | Power across cohort sizes, with `plot()` |
| `cqtc_sample_size()` | Smallest cohort size for a target power (analytical or simulation) |
| `cqtc_simulate()` | Empirical power, coverage and bias for all six models |
| `simulate_sad_data()` | Simulate one study |
| `fit_cqtc_model()` | Fit one model; estimate and confidence bounds at C<sub>supra</sub> |
| `conc_moments()` | Concentration summaries used by the analytical formulas |
| `cqtc_app()` | Launch the Shiny app |

See `?cqtc_power` etc. for details, and the vignette (`vignettes/cqtcpower.Rmd`)
for a worked example.

---

## Methods and assumptions

**Analytical covariance.** The covariance of (θ̂<sub>Z</sub>, θ̂<sub>C</sub>) is derived from
the generalised least-squares information with known variance components. For the
random-slope models, concentrations enter as C*<sub>ik</sub> = C<sub>ik</sub> / (1 + g<sub>i</sub>),
with g<sub>i</sub> = trace(B<sub>i</sub>Σ<sup>−1</sup>), plus a term in
g<sub>i</sub>²/(1 + g<sub>i</sub>)². Concentration summaries are averaged over simulated studies
when concentrations are random.

**Assumptions of the closed form.**

* Variance components are known (no allowance for estimating them).
* The treatment indicator and concentrations are uncorrelated with baseline QTc.
* For `RS0` / `RS0 + Z`, the fixed-effect block of the information matrix is
  approximated by its random-intercept form, and the random intercept and slope
  are uncorrelated.
* `RI` / `RI + Z` ignore any random slope in the data.

In the settings of the accompanying paper, the closed-form SE(Û) is within about
3% of the exact GLS value, and analytical power agrees with simulation within
Monte Carlo error. Because it ignores the uncertainty in the variance components,
analytical power can be slightly optimistic for small studies; check borderline
designs by simulation, ideally with Kenward-Roger degrees of freedom.

**Model fitting.** Simulations fit the models by REML with `lmerTest::lmer()`.
Concentration enters as C/C<sub>supra</sub>, an exact reparameterisation that keeps the
optimisation well scaled; if `bobyqa` reports a problem, `Nelder_Mead` and then
`nloptwrap` are tried. Satterthwaite and Kenward-Roger change only the degrees of
freedom of the t quantile; the standard error is the model-based one.

---

## License

MIT. See [LICENSE](LICENSE).
