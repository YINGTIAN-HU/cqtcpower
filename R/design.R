# ---------------------------------------------------------------------------
# Concentration models, study design and model parameters
# ---------------------------------------------------------------------------

#' Concentrations from a one-compartment PK model
#'
#' Describes plasma concentrations generated from a one-compartment model with
#' first-order absorption, one cohort per dose level. Between-subject
#' variability enters as log-normal random effects on `Ka`, `CL` and `V`.
#'
#' @param dose Numeric vector of dose levels (mg), one per cohort.
#' @param Ka Absorption rate constant (1/h).
#' @param CL Clearance.
#' @param V Volume of distribution. With the defaults, concentrations are in
#'   mg/mL.
#' @param iiv_sd Standard deviation of the log-normal between-subject random
#'   effects on `Ka`, `CL` and `V`. Set to `0` for deterministic profiles.
#' @param resid_cv Standard deviation of an optional multiplicative log-normal
#'   error on each concentration. `0` (default) means no error.
#'
#' @return An object of class `cqtc_conc`.
#' @seealso [conc_profile()], [cqtc_design()]
#' @examples
#' conc_pk_1cmt()
#' conc_pk_1cmt(dose = c(10, 50, 200), iiv_sd = 0.2)
#' @export
conc_pk_1cmt <- function(dose = c(1, 3, 10, 30, 100, 200, 400, 800),
                         Ka = 0.8, CL = 12000, V = 40000,
                         iiv_sd = 0.1, resid_cv = 0) {
  .check_positive(dose, "dose")
  .check_positive(c(Ka, CL, V), "Ka, CL and V")
  .check_nonneg(c(iiv_sd, resid_cv), "iiv_sd and resid_cv")
  if (isTRUE(all.equal(Ka, CL / V))) {
    stop("`Ka` must differ from the elimination rate `CL / V`.", call. = FALSE)
  }
  structure(
    list(type = "pk_1cmt", dose = dose, Ka = Ka, CL = CL, V = V,
         iiv_sd = iiv_sd, resid_cv = resid_cv,
         n_cohort = length(dose),
         random = iiv_sd > 0 || resid_cv > 0),
    class = "cqtc_conc"
  )
}

#' Concentrations from a fixed time profile
#'
#' Describes concentrations as `cmax * cmax_fraction[g] * profile`, i.e. a
#' common time profile (expressed as a fraction of the cohort's maximum
#' concentration) scaled for each cohort `g`. This is the construction used for
#' the sample-size figures of the accompanying paper.
#'
#' @param cmax Maximum mean concentration (reached in the highest cohort).
#' @param cmax_fraction Numeric vector, one entry per cohort: the cohort's
#'   maximum concentration as a fraction of `cmax`.
#' @param profile Numeric vector, one entry per ECG time point: concentration
#'   at that time as a fraction of the cohort's maximum.
#' @param cv Coefficient of variation of normally distributed between-subject
#'   noise added to each concentration. `0` (default) gives identical
#'   concentrations for all active subjects in a cohort.
#'
#' @return An object of class `cqtc_conc`.
#' @seealso [conc_pk_1cmt()], [cqtc_design()]
#' @examples
#' conc_profile()
#' # four cohorts, shorter sampling schedule
#' conc_profile(cmax_fraction = c(1/144, 1/48, 1/16, 1),
#'              profile = c(0.9, 0.95, 0.85, 0.4))
#' @export
conc_profile <- function(cmax = 0.008,
                         cmax_fraction = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1),
                         profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01),
                         cv = 0) {
  .check_positive(cmax, "cmax")
  .check_positive(cmax_fraction, "cmax_fraction")
  .check_nonneg(profile, "profile")
  .check_nonneg(cv, "cv")
  structure(
    list(type = "profile", cmax = cmax, cmax_fraction = cmax_fraction,
         profile = profile, cv = cv,
         n_cohort = length(cmax_fraction),
         random = cv > 0),
    class = "cqtc_conc"
  )
}

#' @export
print.cqtc_conc <- function(x, ...) {
  if (x$type == "pk_1cmt") {
    cat("One-compartment PK concentrations\n")
    cat("  doses (mg):", paste(x$dose, collapse = ", "), "\n")
    cat(sprintf("  Ka = %g, CL = %g, V = %g, IIV sd = %g, residual CV = %g\n",
                x$Ka, x$CL, x$V, x$iiv_sd, x$resid_cv))
  } else {
    cat("Fixed concentration profile\n")
    cat(sprintf("  Cmax = %g; %d cohorts at fractions %s\n", x$cmax, x$n_cohort,
                paste(signif(x$cmax_fraction, 3), collapse = ", ")))
    cat("  time profile:", paste(x$profile, collapse = ", "), "\n")
    if (x$cv > 0) cat("  between-subject CV:", x$cv, "\n")
  }
  invisible(x)
}

#' Parallel SAD study design
#'
#' @param conc Concentration model: [conc_pk_1cmt()] or [conc_profile()]. The
#'   number of cohorts is taken from it.
#' @param n_active,n_placebo Number of active and placebo subjects per cohort.
#'   Their ratio is kept fixed when the sample size is varied, so the cohort
#'   size moves in blocks of `(n_active + n_placebo) / gcd(n_active, n_placebo)`.
#' @param times ECG time points (h). Defaults to `c(1, 2, 4, 8, 12, 16, 24)` for
#'   PK concentrations and to `seq_along(profile)` for a fixed profile, whose
#'   length must match `times`.
#'
#' @return An object of class `cqtc_design`.
#' @examples
#' cqtc_design()
#' cqtc_design(conc_profile(), n_active = 6, n_placebo = 2)
#' @export
cqtc_design <- function(conc = conc_pk_1cmt(), n_active = 6, n_placebo = 2,
                        times = NULL) {
  if (!inherits(conc, "cqtc_conc")) {
    stop("`conc` must be created with conc_pk_1cmt() or conc_profile().",
         call. = FALSE)
  }
  .check_count(n_active, "n_active", min = 1)
  .check_count(n_placebo, "n_placebo", min = 0)
  if (is.null(times)) {
    times <- if (conc$type == "pk_1cmt") c(1, 2, 4, 8, 12, 16, 24) else seq_along(conc$profile)
  }
  .check_nonneg(times, "times")
  if (conc$type == "profile" && length(times) != length(conc$profile)) {
    stop("`times` must have one entry per element of the concentration `profile`.",
         call. = FALSE)
  }
  g <- .gcd(n_active, n_placebo)
  structure(
    list(conc = conc, n_active = as.integer(n_active),
         n_placebo = as.integer(n_placebo), times = times,
         n_cohort = conc$n_cohort, n_time = length(times),
         block = as.integer((n_active + n_placebo) / g)),
    class = "cqtc_design"
  )
}

#' @export
print.cqtc_design <- function(x, ...) {
  cat("Parallel SAD design\n")
  cat(sprintf("  %d cohorts x (%d active + %d placebo) = %d subjects\n",
              x$n_cohort, x$n_active, x$n_placebo,
              x$n_cohort * (x$n_active + x$n_placebo)))
  cat("  ECG times:", paste(x$times, collapse = ", "), "\n")
  cat("  ")
  print(x$conc)
  invisible(x)
}

#' True effects and variance components
#'
#' Collects the quantities that determine power: the true \eqn{\Delta\Delta QTc}
#' at the supratherapeutic concentration, the variance components of the
#' mixed-effects model, and the decision rule.
#'
#' The true concentration slope is `theta_C = (U1 - theta_Z) / C_supra`. The
#' random-slope standard deviation is `sigma1 = slope_cv * |theta_C|` unless
#' `sigma1` is given directly.
#'
#' @param U1 True \eqn{\Delta\Delta QTc} (ms) at `C_supra`:
#'   `theta_Z + theta_C * C_supra`.
#' @param C_supra Supratherapeutic concentration (same units as the
#'   concentrations).
#' @param theta_Z True treatment (drug, concentration-independent) effect (ms).
#' @param sigma_e Residual standard deviation (ms).
#' @param sigma0 Standard deviation of the random intercept (ms).
#' @param slope_cv Coefficient of variation of the random slope relative to
#'   `theta_C`. `0` means no random slope in the data.
#' @param sigma1 Standard deviation of the random slope (ms per concentration
#'   unit). Overrides `slope_cv` when supplied.
#' @param rho Correlation between random intercept and random slope (used by
#'   the simulation only; the analytical formulas assume `rho = 0`).
#' @param U0 Regulatory threshold (ms). "No QT effect" is concluded when the
#'   upper confidence bound is below `U0`.
#' @param alpha One-sided significance level.
#' @param theta0,theta_W,theta_T Intercept, baseline effect and time effects
#'   used when simulating data. They do not affect power. `theta_T` defaults to
#'   the values in the paper for seven time points and to zeros otherwise.
#' @param baseline_mean,baseline_sd Mean and SD of baseline QTc (ms), used when
#'   simulating data.
#'
#' @return An object of class `cqtc_params`.
#' @examples
#' cqtc_params(U1 = 5)
#' cqtc_params(U1 = 5, theta_Z = 0.5, slope_cv = 0.75)
#' cqtc_params(U1 = 3, sigma1 = 600)
#' @export
cqtc_params <- function(U1, C_supra = 0.008, theta_Z = 0,
                        sigma_e = 5, sigma0 = 8,
                        slope_cv = 0, sigma1 = NULL, rho = 0,
                        U0 = 10, alpha = 0.05,
                        theta0 = -5.76, theta_W = -0.58, theta_T = NULL,
                        baseline_mean = 408, baseline_sd = 16) {
  .check_scalar(U1, "U1")
  .check_positive(C_supra, "C_supra")
  .check_scalar(theta_Z, "theta_Z")
  .check_positive(c(sigma_e), "sigma_e")
  .check_nonneg(c(sigma0, slope_cv), "sigma0 and slope_cv")
  if (!is.numeric(rho) || length(rho) != 1 || abs(rho) > 1) {
    stop("`rho` must be a single number in [-1, 1].", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 0.5) {
    stop("`alpha` must be a single number in (0, 0.5).", call. = FALSE)
  }
  theta_C <- (U1 - theta_Z) / C_supra
  if (is.null(sigma1)) {
    sigma1 <- slope_cv * abs(theta_C)
  } else {
    .check_nonneg(sigma1, "sigma1")
    slope_cv <- if (theta_C != 0) sigma1 / abs(theta_C) else NA_real_
  }
  structure(
    list(U1 = U1, C_supra = C_supra, theta_Z = theta_Z, theta_C = theta_C,
         sigma_e = sigma_e, sigma0 = sigma0, sigma1 = sigma1,
         slope_cv = slope_cv, rho = rho, U0 = U0, alpha = alpha,
         theta0 = theta0, theta_W = theta_W, theta_T = theta_T,
         baseline_mean = baseline_mean, baseline_sd = baseline_sd),
    class = "cqtc_params"
  )
}

#' @export
print.cqtc_params <- function(x, ...) {
  cat("C-QTc parameters\n")
  cat(sprintf("  true ddQTc at C_supra = %g: U1 = %g ms (theta_Z = %g, theta_C = %.4g)\n",
              x$C_supra, x$U1, x$theta_Z, x$theta_C))
  cat(sprintf("  sigma_e = %g, sigma0 = %g, sigma1 = %.4g (CV %s), rho = %g\n",
              x$sigma_e, x$sigma0, x$sigma1,
              if (is.na(x$slope_cv)) "NA" else sprintf("%.0f%%", 100 * x$slope_cv),
              x$rho))
  cat(sprintf("  decision: one-sided %.0f%% upper bound < U0 = %g ms\n",
              100 * (1 - x$alpha), x$U0))
  invisible(x)
}

# ---- allocation helpers -----------------------------------------------------

# Translate a cohort size into (n_active, n_placebo) keeping the design ratio.
.allocation <- function(design, n_per_cohort = NULL) {
  if (is.null(n_per_cohort)) {
    return(list(n_active = design$n_active, n_placebo = design$n_placebo))
  }
  .check_count(n_per_cohort, "n_per_cohort", min = 1)
  if (n_per_cohort %% design$block != 0) {
    stop(sprintf(
      "`n_per_cohort` must be a multiple of %d to keep the %d:%d active:placebo ratio.",
      design$block, design$n_active, design$n_placebo), call. = FALSE)
  }
  k <- n_per_cohort / (design$n_active + design$n_placebo)
  list(n_active = as.integer(round(design$n_active * k)),
       n_placebo = as.integer(round(design$n_placebo * k)))
}

.theta_T <- function(params, n_time) {
  th <- params$theta_T
  if (is.null(th)) {
    th <- if (n_time == 7) c(1.93, 2.29, 3.53, -3.15, -0.28, 5.17, 3.98) else rep(0, n_time)
  }
  if (length(th) != n_time) {
    stop("`theta_T` must have one entry per ECG time point.", call. = FALSE)
  }
  th
}

# ---- argument checks ----------------------------------------------------------

.check_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
    stop(sprintf("`%s` must be a single finite number.", name), call. = FALSE)
  }
}
.check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0 || any(!is.finite(x)) || any(x <= 0)) {
    stop(sprintf("`%s` must be positive.", name), call. = FALSE)
  }
}
.check_nonneg <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0 || any(!is.finite(x)) || any(x < 0)) {
    stop(sprintf("`%s` must be non-negative.", name), call. = FALSE)
  }
}
.check_count <- function(x, name, min = 0) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x != round(x) || x < min) {
    stop(sprintf("`%s` must be a whole number >= %d.", name, min), call. = FALSE)
  }
}
.gcd <- function(a, b) if (b == 0) a else .gcd(b, a %% b)
