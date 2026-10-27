#' cqtcpower: power and sample size for C-QTc analyses in parallel SAD studies
#'
#' The package answers two questions for a planned parallel-group single
#' ascending dose (SAD) study whose ECG data will be analysed with a
#' concentration-QTc (C-QTc) linear mixed-effects model:
#'
#' * **What is the power** to conclude "no QT effect" (upper bound of the
#'   one-sided \eqn{1-\alpha} confidence interval of the predicted
#'   \eqn{\Delta\Delta QTc} at the supratherapeutic concentration below
#'   \eqn{U_0 = 10} ms) for a given number of subjects per cohort?
#'   See [cqtc_power()] and [cqtc_power_curve()].
#' * **How many subjects per cohort** are needed to reach a target power?
#'   See [cqtc_sample_size()].
#'
#' Analytical results use the closed-form covariance approximations for the
#' `RI`, `RI + Z`, `RS0` and `RS0 + Z` models. Simulation results
#' ([cqtc_simulate()]) are available for all six models, including the
#' correlated random-slope models `RS` and `RS + Z`.
#'
#' A typical workflow:
#' 1. describe the concentrations with [conc_pk_1cmt()] or [conc_profile()];
#' 2. describe the study with [cqtc_design()];
#' 3. set the true effect and variance components with [cqtc_params()];
#' 4. call [cqtc_power()], [cqtc_power_curve()], [cqtc_sample_size()] or
#'    [cqtc_simulate()].
#'
#' @section Models:
#' Every model contains fixed effects for baseline QTc and for each time point
#' (no global intercept).
#'
#' | Model     | Fixed effects              | Random effects (by subject)          |
#' |-----------|----------------------------|--------------------------------------|
#' | `RI`      | concentration              | intercept                            |
#' | `RI + Z`  | treatment + concentration  | intercept                            |
#' | `RS0`     | concentration              | intercept + slope, uncorrelated      |
#' | `RS0 + Z` | treatment + concentration  | intercept + slope, uncorrelated      |
#' | `RS`      | concentration              | intercept + slope, correlated        |
#' | `RS + Z`  | treatment + concentration  | intercept + slope, correlated        |
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats pnorm qnorm qt rnorm var cov setNames sd
NULL
