# ---------------------------------------------------------------------------
# Analytical covariance, power and sample size (Appendix A of the paper)
# ---------------------------------------------------------------------------

.analytic_models <- c("RI", "RI + Z", "RS0", "RS0 + Z")
.all_models <- c("RI", "RI + Z", "RS0", "RS0 + Z", "RS", "RS + Z")

# Accept "RI+Z", "ri + z", ... and return canonical names.
.match_models <- function(models, allowed = .all_models) {
  key <- function(x) toupper(gsub("\\s+", "", x))
  idx <- match(key(models), key(allowed))
  if (anyNA(idx)) {
    stop(sprintf("Unknown model(s): %s. Choose from: %s.",
                 paste(models[is.na(idx)], collapse = ", "),
                 paste(allowed, collapse = ", ")), call. = FALSE)
  }
  allowed[idx]
}
.has_Z <- function(model) grepl("Z", model, fixed = TRUE)
.has_slope <- function(model) grepl("^RS", model)

# g_i = trace(B_i Sigma^{-1}) for each subject (columns of C).
.gi <- function(C, sigma_e, sigma0, sigma1) {
  T <- nrow(C)
  u <- sigma_e^2 + T * sigma0^2
  sigma1^2 / sigma_e^2 * colSums(C^2) -
    sigma0^2 * sigma1^2 / (sigma_e^2 * u) * colSums(C)^2
}

# Concentration summaries entering the covariance formulas, for one study.
# C* = C / (1 + g_i), as in Appendix A.2.
.conc_moments_one <- function(C, sigma_e, sigma0, sigma1) {
  g <- .gi(C, sigma_e, sigma0, sigma1)
  Cs <- sweep(C, 2, 1 + g, "/")
  list(mean_C = mean(colSums(C)), var_C = var(colSums(C)),
       var_CT = apply(C, 1, var),
       mean_Cs = mean(colSums(Cs)), var_Cs = var(colSums(Cs)),
       var_CTs = apply(Cs, 1, var),
       mean_gi2 = mean(g^2 / (1 + g)^2))
}

#' Concentration summaries used by the analytical formulas
#'
#' Averages, over `n_sim` simulated studies, the concentration moments that
#' enter the closed-form covariance of \eqn{(\hat\theta_Z, \hat\theta_C)}:
#' the mean and variance of subject sums, the variance at each time point, and
#' the same quantities for \eqn{C^*_{ik} = C_{ik} / (1 + g_i)}, plus
#' \eqn{\overline{g_i^2/(1+g_i)^2}}. Deterministic concentrations need only one
#' study.
#'
#' @inheritParams cqtc_power
#' @return A named list of moments.
#' @keywords internal
#' @export
conc_moments <- function(design, params, n_conc_sim = 1000, seed = 1) {
  stopifnot(inherits(design, "cqtc_design"), inherits(params, "cqtc_params"))
  n_sim <- if (design$conc$random) n_conc_sim else 1
  .with_seed(seed, {
    res <- lapply(seq_len(n_sim), function(i) {
      st <- .study_conc(design, design$n_active, design$n_placebo)
      .conc_moments_one(st$C, params$sigma_e, params$sigma0, params$sigma1)
    })
  })
  avg <- function(name) {
    v <- lapply(res, `[[`, name)
    if (length(v[[1]]) > 1) rowMeans(do.call(cbind, v)) else mean(unlist(v))
  }
  stats::setNames(lapply(names(res[[1]]), avg), names(res[[1]]))
}

# Inverse covariance (information) of (theta_Z, theta_C) for N subjects.
# Random-slope models use C* and the g^2/(1+g)^2 term; random-intercept models
# use C (sigma1 = 0).
.info_matrix <- function(mom, N, R, T, sigma_e, sigma0, sigma1, slope) {
  u <- sigma_e^2 + T * sigma0^2
  slope <- slope && sigma1 > 0
  mC  <- if (slope) mom$mean_Cs else mom$mean_C
  vC  <- if (slope) mom$var_Cs  else mom$var_C
  vCT <- if (slope) mom$var_CTs else mom$var_CT
  M11 <- T / u * (R^2 + R) / (R + 1)^3 * N
  M12 <- 1 / (R + 1) * mC / u * N
  M22 <- (N - 1) * sum(vCT) / sigma_e^2 - sigma0^2 / (sigma_e^2 * u) * vC * (N - 1)
  if (slope) M22 <- M22 + mom$mean_gi2 * N / sigma1^2
  matrix(c(M11, M12, M12, M22), 2, 2)
}

# Standard error of U = theta_Z + theta_C * C_supra (or theta_C * C_supra).
.se_U <- function(M, C_supra, with_Z) {
  if (with_Z) {
    cc <- c(1, C_supra)
    sqrt(drop(crossprod(cc, solve(M, cc))))
  } else {
    C_supra / sqrt(M[2, 2])
  }
}

.power_from_se <- function(se, params) {
  pnorm((params$U0 - params$U1) / se - qnorm(1 - params$alpha))
}

.analytic_table <- function(design, params, models, n_grid, mom) {
  R <- design$n_active / design$n_placebo
  T <- design$n_time
  rows <- lapply(n_grid, function(n) {
    al <- .allocation(design, n)
    N <- design$n_cohort * (al$n_active + al$n_placebo)
    lapply(models, function(m) {
      M <- .info_matrix(mom, N, R, T, params$sigma_e, params$sigma0,
                        params$sigma1, slope = .has_slope(m))
      se <- .se_U(M, params$C_supra, .has_Z(m))
      data.frame(model = m, n_active = al$n_active, n_placebo = al$n_placebo,
                 n_per_cohort = al$n_active + al$n_placebo, n_total = N,
                 se_U = se, power = .power_from_se(se, params),
                 stringsAsFactors = FALSE)
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}

.check_analytic_inputs <- function(design, params, models) {
  if (!inherits(design, "cqtc_design")) stop("`design` must come from cqtc_design().", call. = FALSE)
  if (!inherits(params, "cqtc_params")) stop("`params` must come from cqtc_params().", call. = FALSE)
  models <- .match_models(models)
  bad <- setdiff(models, .analytic_models)
  if (length(bad)) {
    stop(sprintf(paste0("No analytical formula for %s (correlated random slope). ",
                        "Use cqtc_simulate(), or cqtc_sample_size(method = \"simulation\")."),
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  if (design$n_placebo == 0) {
    stop("The analytical formulas need placebo subjects (n_placebo >= 1).", call. = FALSE)
  }
  models
}

#' Analytical power of a C-QTc analysis
#'
#' Power to conclude "no QT effect": the upper bound of the one-sided
#' \eqn{1-\alpha} confidence interval of
#' \eqn{U = \theta_Z + \theta_C C_{supra}} (or \eqn{\theta_C C_{supra}} for
#' models without a treatment effect) falls below `U0`:
#' \deqn{\mathrm{power} = \Phi\left(\frac{U_0 - U_1}{\mathrm{SE}(\hat U)} - z_{1-\alpha}\right),}
#' with \eqn{\mathrm{SE}(\hat U)} from the closed-form covariance
#' approximations (assuming the treatment indicator and the concentrations are
#' uncorrelated with baseline QTc).
#'
#' `RI` and `RI + Z` use the random-intercept covariance (ignoring any random
#' slope in the data); `RS0` and `RS0 + Z` use the uncorrelated random-slope
#' covariance with `sigma1` from `params`. For the correlated random-slope
#' models `RS` and `RS + Z` use [cqtc_simulate()].
#'
#' For random PK concentrations the concentration moments are averaged over
#' `n_conc_sim` simulated studies (no model fitting, so this is fast).
#'
#' @param design A [cqtc_design()].
#' @param params A [cqtc_params()].
#' @param models Models to evaluate; any of `"RI"`, `"RI + Z"`, `"RS0"`,
#'   `"RS0 + Z"`.
#' @param n_per_cohort Optional vector of cohort sizes (subjects per cohort) at
#'   which to evaluate power, keeping the active:placebo ratio of `design`.
#'   Defaults to the design's cohort size.
#' @param n_conc_sim Number of simulated studies used to average the
#'   concentration moments when concentrations are random.
#' @param seed Random seed for the concentration simulation (`NULL` to use the
#'   current RNG state). The global RNG state is restored afterwards.
#'
#' @return A data frame of class `cqtc_power` with one row per model and cohort
#'   size: allocation, total sample size, `se_U` and `power`.
#' @seealso [cqtc_power_curve()], [cqtc_sample_size()], [cqtc_simulate()]
#' @examples
#' des <- cqtc_design(conc_pk_1cmt(), n_active = 6, n_placebo = 2)
#' cqtc_power(des, cqtc_params(U1 = 5, theta_Z = 0.5, slope_cv = 0.75))
#' @export
cqtc_power <- function(design, params, models = .analytic_models,
                       n_per_cohort = NULL, n_conc_sim = 1000, seed = 1) {
  models <- .check_analytic_inputs(design, params, models)
  if (is.null(n_per_cohort)) n_per_cohort <- design$n_active + design$n_placebo
  for (n in n_per_cohort) .allocation(design, n)
  mom <- conc_moments(design, params, n_conc_sim = n_conc_sim, seed = seed)
  tab <- .analytic_table(design, params, models, n_per_cohort, mom)
  tab <- tab[order(match(tab$model, models), tab$n_per_cohort), ]
  rownames(tab) <- NULL
  structure(tab, class = c("cqtc_power", "data.frame"),
            design = design, params = params)
}

#' Power as a function of cohort size
#'
#' Evaluates [cqtc_power()] over a grid of cohort sizes. The result has a
#' `plot()` method.
#'
#' @inheritParams cqtc_power
#' @param n_per_cohort Cohort sizes to evaluate. Defaults to 1 to 10 blocks of
#'   the design's active:placebo ratio (e.g. 4, 8, ..., 40 for 6:2).
#' @param target Optional target power drawn as a reference line by `plot()`.
#' @return A `cqtc_power` data frame (see [cqtc_power()]).
#' @examples
#' pc <- cqtc_power_curve(cqtc_design(conc_profile()),
#'                        cqtc_params(U1 = 5, sigma1 = 600))
#' head(pc)
#' plot(pc)
#' @export
cqtc_power_curve <- function(design, params, models = .analytic_models,
                             n_per_cohort = NULL, target = 0.8,
                             n_conc_sim = 1000, seed = 1) {
  if (is.null(n_per_cohort)) n_per_cohort <- design$block * seq_len(10)
  out <- cqtc_power(design, params, models = models, n_per_cohort = n_per_cohort,
                    n_conc_sim = n_conc_sim, seed = seed)
  attr(out, "target") <- target
  out
}

#' Sample size for a C-QTc analysis
#'
#' Smallest cohort size (keeping the active:placebo ratio of `design`) for
#' which the power reaches `power`.
#'
#' With `method = "analytical"` (the default) the closed-form power of
#' [cqtc_power()] is used. A first guess comes from
#' \deqn{N = \left(\frac{z_{1-\beta} + z_{1-\alpha}}{U_0 - U_1}\right)^2 N_{ref}\,\mathrm{SE}^2_{ref}(\hat U),}
#' rounded up to whole blocks per cohort; the search then moves block by block
#' to the smallest cohort size whose analytical power is at least `power`.
#' `n_formula` reports the rounded closed-form value.
#'
#' With `method = "simulation"`, empirical power from [cqtc_simulate()] is used
#' instead, starting from the analytical answer (for `RS` and `RS + Z`, from
#' `RS0` and `RS0 + Z`). This is the only option for the correlated
#' random-slope models and can take a long time: every step fits `nsim`
#' mixed models.
#'
#' @inheritParams cqtc_power
#' @param power Target power.
#' @param models Models to size for. With `method = "simulation"`, `"RS"` and
#'   `"RS + Z"` are also allowed.
#' @param method `"analytical"` or `"simulation"`.
#' @param max_n_per_cohort Upper limit of the search.
#' @param nsim,ddf,parallel Passed to [cqtc_simulate()] when
#'   `method = "simulation"`. `ddf` must be a single method here.
#'
#' @return A data frame of class `cqtc_sample_size` with, per model, the
#'   allocation per cohort, total sample size, the power achieved and its
#'   standard error, and `n_formula` (analytical method).
#' @examples
#' des <- cqtc_design(conc_profile(), n_active = 6, n_placebo = 2)
#' cqtc_sample_size(des, cqtc_params(U1 = 5, sigma1 = 600), power = 0.8)
#' @export
cqtc_sample_size <- function(design, params, power = 0.8,
                             models = .analytic_models,
                             method = c("analytical", "simulation"),
                             max_n_per_cohort = 500,
                             n_conc_sim = 1000, seed = 1,
                             nsim = 500, ddf = "Kenward-Roger", parallel = FALSE) {
  method <- match.arg(method)
  if (!is.numeric(power) || length(power) != 1 || power <= 0 || power >= 1) {
    stop("`power` must be a single number in (0, 1).", call. = FALSE)
  }
  if (params$U1 >= params$U0) {
    stop("The true ddQTc `U1` must be below the threshold `U0`; power cannot reach the target.",
         call. = FALSE)
  }
  if (method == "simulation") {
    return(.sample_size_sim(design, params, power, .match_models(models),
                            max_n_per_cohort, n_conc_sim, seed, nsim, ddf, parallel))
  }
  models <- .check_analytic_inputs(design, params, models)
  mom <- conc_moments(design, params, n_conc_sim = n_conc_sim, seed = seed)
  b <- design$block
  n_ref <- design$n_active + design$n_placebo
  z <- qnorm(power) + qnorm(1 - params$alpha)
  delta <- params$U0 - params$U1
  rows <- lapply(models, function(m) {
    ref <- .analytic_table(design, params, m, n_ref, mom)
    N_req <- (z / delta)^2 * ref$se_U^2 * ref$n_total
    n_formula <- max(b, ceiling(N_req / (design$n_cohort * b)) * b)
    pw <- function(n) .analytic_table(design, params, m, n, mom)$power
    n <- min(n_formula, max_n_per_cohort - max_n_per_cohort %% b)
    while (pw(n) < power && n + b <= max_n_per_cohort) n <- n + b
    while (n - b >= b && pw(n - b) >= power) n <- n - b
    res <- .analytic_table(design, params, m, n, mom)
    if (res$power < power) {
      warning(sprintf("%s: target power not reached within max_n_per_cohort = %d.",
                      m, max_n_per_cohort), call. = FALSE)
    }
    res$n_formula <- n_formula
    res
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("cqtc_sample_size", "data.frame"),
            design = design, params = params, target = power, method = method)
}
