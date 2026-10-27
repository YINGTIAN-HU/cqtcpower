# ---------------------------------------------------------------------------
# Data simulation, model fitting and empirical power
# ---------------------------------------------------------------------------

#' Simulate one parallel SAD study
#'
#' Generates concentrations and \eqn{\Delta QTc} for every subject and time
#' point from
#' \deqn{\Delta QTc_{ik} = \theta_0 + \eta_{0i} + \theta_{Tk} + \theta_W W_i +
#'   \theta_Z Z_i + (\theta_C + \eta_{1i}) C_{ik} + \epsilon_{ik},}
#' with \eqn{(\eta_{0i}, \eta_{1i})} bivariate normal (SDs `sigma0`, `sigma1`,
#' correlation `rho`), \eqn{\epsilon_{ik} \sim N(0, \sigma_e^2)}, \eqn{W_i} the
#' centred baseline QTc and \eqn{Z_i} the active-treatment indicator.
#'
#' @inheritParams cqtc_power
#' @param n_per_cohort Optional cohort size (keeps the design's ratio).
#' @param seed Optional random seed (the global RNG state is restored).
#' @return A data frame in long format with columns `Subject`, `Cohort`,
#'   `Treatment`, `Time` (factor), `Concentration`, `Baseline` (centred) and
#'   `delta_QTC`; plus `Dose` for PK designs.
#' @examples
#' d <- simulate_sad_data(cqtc_design(), cqtc_params(U1 = 5), seed = 1)
#' head(d)
#' @export
simulate_sad_data <- function(design, params, n_per_cohort = NULL, seed = NULL) {
  stopifnot(inherits(design, "cqtc_design"), inherits(params, "cqtc_params"))
  al <- .allocation(design, n_per_cohort)
  .with_seed(seed, .simulate_sad(design, params, al$n_active, al$n_placebo))
}

.simulate_sad <- function(design, params, n_active, n_placebo) {
  T <- design$n_time
  theta_T <- .theta_T(params, T)
  st <- .study_conc(design, n_active, n_placebo)
  N <- ncol(st$C)
  W <- rnorm(N, params$baseline_mean, params$baseline_sd) - params$baseline_mean
  s0 <- params$sigma0
  s1 <- params$sigma1
  Sigma <- matrix(c(s0^2, params$rho * s0 * s1, params$rho * s0 * s1, s1^2), 2)
  eta <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma)
  eta <- matrix(eta, ncol = 2)
  conc <- as.vector(st$C)                       # subject-major (time within subject)
  subj <- rep(seq_len(N), each = T)
  y <- params$theta0 + rep(theta_T, N) + params$theta_W * W[subj] +
    params$theta_Z * st$treatment[subj] +
    (params$theta_C + eta[subj, 2]) * conc + eta[subj, 1] +
    rnorm(N * T, 0, params$sigma_e)
  d <- data.frame(Subject = factor(subj), Cohort = st$cohort[subj],
                  Treatment = st$treatment[subj],
                  Time = factor(rep(design$times, N), levels = design$times),
                  Concentration = conc, Baseline = W[subj], delta_QTC = y)
  if (design$conc$type == "pk_1cmt") {
    d$Dose <- design$conc$dose[d$Cohort] * d$Treatment
  }
  d
}

.ddf_names <- c(Wald = "Wald", Satterthwaite = "Satterthwaite",
                `Kenward-Roger` = "Kenward-Roger")
.match_ddf <- function(ddf) {
  key <- tolower(gsub("[^A-Za-z]", "", ddf))
  map <- c(wald = "Wald", normal = "Wald", z = "Wald",
           satterthwaite = "Satterthwaite", sw = "Satterthwaite",
           kenwardroger = "Kenward-Roger", kr = "Kenward-Roger")
  out <- unname(map[key])
  if (anyNA(out)) {
    stop("`ddf` must be any of \"Wald\", \"Satterthwaite\", \"Kenward-Roger\".", call. = FALSE)
  }
  if ("Kenward-Roger" %in% out && !requireNamespace("pbkrtest", quietly = TRUE)) {
    warning("Package 'pbkrtest' is not installed; dropping Kenward-Roger.", call. = FALSE)
    out <- setdiff(out, "Kenward-Roger")
    if (!length(out)) out <- "Satterthwaite"
  }
  unique(out)
}

.model_formula <- function(model) {
  # C_rel = Concentration / C_supra: the slope is the ddQTc at C_supra and the
  # random-slope SD is on the ms scale, which keeps the optimisation well scaled.
  fixed <- if (.has_Z(model)) {
    "delta_QTC ~ 0 + Time + Baseline + Treatment + C_rel"
  } else {
    "delta_QTC ~ 0 + Time + Baseline + C_rel"
  }
  random <- switch(sub(" \\+ Z$", "", model),
                   RI  = "(1 | Subject)",
                   RS0 = "(C_rel || Subject)",
                   RS  = "(C_rel | Subject)")
  stats::as.formula(paste(fixed, "+", random))
}

.optimizers <- c("bobyqa", "Nelder_Mead", "nloptwrap")

# Fit with one optimizer; collect (do not print) warnings.
.fit_lmer <- function(formula, data, optimizer) {
  msgs <- character(0)
  fit <- withCallingHandlers(
    suppressMessages(lmerTest::lmer(formula, data = data, REML = TRUE,
                                    control = lme4::lmerControl(optimizer = optimizer))),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, converged = length(msgs) == 0, warnings = msgs, optimizer = optimizer)
}

#' Fit a C-QTc mixed model and estimate ddQTc at C_supra
#'
#' Fits one of the six linear mixed-effects models with `lmerTest::lmer()`
#' (REML) and returns the estimate of
#' \eqn{U = \theta_Z + \theta_C C_{supra}} (or \eqn{\theta_C C_{supra}}), its
#' standard error and two-sided \eqn{1 - 2\alpha} confidence limits (whose upper
#' limit is the one-sided \eqn{1-\alpha} bound) for each degrees-of-freedom
#' method. The standard error is the model-based one; Satterthwaite and
#' Kenward-Roger only change the degrees of freedom of the t quantile.
#'
#' **Convergence.** Concentration enters the model as
#' `C_rel = Concentration / C_supra`. This is an exact reparameterisation (the
#' fitted model and \eqn{\hat U} are unchanged) but puts the random-slope SD on
#' the same scale as the other variance components, which avoids most
#' convergence failures of random-slope fits on the mg/mL scale. If `bobyqa`
#' still reports a problem, the fit is repeated with `Nelder_Mead` and then
#' `nloptwrap`; the first clean fit is kept, otherwise the one with the highest
#' REML log-likelihood. Boundary (singular) fits, where a variance component is
#' estimated as zero, are valid estimates and are flagged, not refitted.
#'
#' @param data Data in the format of [simulate_sad_data()]: columns
#'   `delta_QTC`, `Time` (factor), `Baseline`, `Concentration`, `Subject` and,
#'   for models with Z, `Treatment`.
#' @param model One of `"RI"`, `"RI + Z"`, `"RS0"`, `"RS0 + Z"`, `"RS"`,
#'   `"RS + Z"`.
#' @param C_supra Supratherapeutic concentration.
#' @param alpha One-sided significance level.
#' @param ddf Degrees-of-freedom methods: any of `"Wald"` (normal quantile),
#'   `"Satterthwaite"`, `"Kenward-Roger"` (needs 'pbkrtest').
#' @param keep_fit Return the fitted model as attribute `"fit"`.
#'
#' @return A one-row data frame with `estimate`, `se`, `converged` (no
#'   optimiser warning), `singular` (boundary fit), `optimizer` and, per method,
#'   `lower_<method>`, `upper_<method>` and `df_<method>`. The optimiser
#'   warnings, if any, are in attribute `"warnings"`.
#' @examples
#' d <- simulate_sad_data(cqtc_design(), cqtc_params(U1 = 5), seed = 1)
#' fit_cqtc_model(d, "RI + Z", C_supra = 0.008, ddf = c("Wald", "Satterthwaite"))
#' @export
fit_cqtc_model <- function(data, model, C_supra, alpha = 0.05,
                           ddf = c("Wald", "Satterthwaite", "Kenward-Roger"),
                           keep_fit = FALSE) {
  model <- .match_models(model)
  if (length(model) != 1) stop("`model` must be a single model.", call. = FALSE)
  ddf <- .match_ddf(ddf)
  data$Time <- factor(data$Time)
  data$C_rel <- data$Concentration / C_supra
  form <- .model_formula(model)
  best <- NULL
  for (opt in .optimizers) {
    r <- .fit_lmer(form, data, opt)
    if (is.null(best) || (!best$converged && stats::logLik(r$fit) > stats::logLik(best$fit))) best <- r
    if (r$converged) { best <- r; break }
  }
  fit <- best$fit
  beta <- lme4::fixef(fit)
  L <- stats::setNames(numeric(length(beta)), names(beta))
  L["C_rel"] <- 1
  if (.has_Z(model)) L["Treatment"] <- 1
  est <- sum(L * beta)
  se <- sqrt(drop(crossprod(L, as.matrix(stats::vcov(fit)) %*% L)))
  out <- data.frame(model = model, estimate = est, se = se,
                    converged = best$converged, singular = lme4::isSingular(fit),
                    optimizer = best$optimizer, stringsAsFactors = FALSE)
  for (m in ddf) {
    df <- switch(m,
      Wald = Inf,
      Satterthwaite = lmerTest::contest(fit, L, joint = FALSE, ddf = "Satterthwaite")$df,
      `Kenward-Roger` = lmerTest::contest(fit, L, joint = FALSE, ddf = "Kenward-Roger")$df)
    q <- if (is.finite(df)) qt(1 - alpha, df) else qnorm(1 - alpha)
    tag <- .ddf_tag(m)
    out[[paste0("df_", tag)]] <- df
    out[[paste0("lower_", tag)]] <- est - q * se
    out[[paste0("upper_", tag)]] <- est + q * se
  }
  if (keep_fit) attr(out, "fit") <- fit
  if (length(best$warnings)) attr(out, "warnings") <- best$warnings
  out
}
.ddf_tag <- function(m) c(Wald = "Wald", Satterthwaite = "SW", `Kenward-Roger` = "KR")[[m]]

#' Simulation-based power, coverage and bias
#'
#' Simulates `nsim` studies with [simulate_sad_data()], fits each requested
#' model to every simulated study with [fit_cqtc_model()] (common random
#' numbers across models), and summarises:
#'
#' * `power`: proportion of fits whose one-sided upper bound is below `U0`;
#' * `coverage`: proportion of two-sided \eqn{1-2\alpha} intervals containing
#'   the true `U1`;
#' * `rel_bias`: mean of \eqn{(\hat U - U_1)/U_1} (`NA` when `U1 = 0`), and
#'   `bias`, the mean of \eqn{\hat U - U_1};
#' * `power_mcse`: Monte Carlo standard error of `power`;
#' * `analytical`: the analytical power from [cqtc_power()] where available.
#'
#' Replicate `i` uses seed `seed + i - 1`, so results do not depend on whether
#' the simulation runs in parallel. Fits that fail with an error are dropped and
#' counted in `n_failed`. `n_nonconv` counts fits for which every optimiser
#' still reported a convergence problem, and `n_singular` boundary fits (see
#' [fit_cqtc_model()]); both are kept unless `drop_nonconverged = TRUE`.
#'
#' @inheritParams cqtc_power
#' @param nsim Number of simulated studies.
#' @param models Models to fit (default: all six).
#' @param ddf Degrees-of-freedom methods for the confidence bounds; see
#'   [fit_cqtc_model()].
#' @param n_per_cohort Optional cohort size (keeps the design's ratio).
#' @param seed Base seed.
#' @param parallel Use `future.apply::future_lapply()` (set up workers with
#'   `future::plan()` first).
#' @param keep_replicates Keep the per-replicate estimates in the result.
#' @param progress Print a dot every 10% of replicates (sequential runs only).
#' @param drop_nonconverged Exclude fits that did not converge from the
#'   summaries (they are counted in `n_nonconv` either way).
#'
#' @return An object of class `cqtc_sim`: a list with `summary` (one row per
#'   model and ddf method), `replicates` (if kept), `design`, `params`, `nsim`.
#' @examples
#' \donttest{
#' sim <- cqtc_simulate(cqtc_design(), cqtc_params(U1 = 5, theta_Z = 0.5),
#'                      nsim = 20, models = c("RI", "RI + Z"), ddf = "Wald")
#' sim
#' }
#' @export
cqtc_simulate <- function(design, params, nsim = 1000, models = .all_models,
                          ddf = c("Wald", "Satterthwaite", "Kenward-Roger"),
                          n_per_cohort = NULL, seed = 1, parallel = FALSE,
                          keep_replicates = TRUE, progress = FALSE,
                          drop_nonconverged = FALSE, n_conc_sim = 1000) {
  stopifnot(inherits(design, "cqtc_design"), inherits(params, "cqtc_params"))
  .check_count(nsim, "nsim", min = 1)
  models <- .match_models(models)
  ddf <- .match_ddf(ddf)
  al <- .allocation(design, n_per_cohort)
  one_rep <- function(i) {
    .with_seed(seed + i - 1, {
      d <- .simulate_sad(design, params, al$n_active, al$n_placebo)
    })
    res <- lapply(models, function(m) {
      tryCatch(fit_cqtc_model(d, m, params$C_supra, params$alpha, ddf),
               error = function(e) NULL)
    })
    res <- do.call(rbind, res)
    if (!is.null(res)) res$rep <- i
    res
  }
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("`parallel = TRUE` needs the 'future.apply' package.", call. = FALSE)
    }
    reps <- future.apply::future_lapply(seq_len(nsim), one_rep, future.seed = TRUE)
  } else {
    step <- max(1, nsim %/% 10)
    reps <- lapply(seq_len(nsim), function(i) {
      r <- one_rep(i)
      if (progress && i %% step == 0) cat(".")
      r
    })
    if (progress) cat("\n")
  }
  reps <- do.call(rbind, reps)
  summary <- .summarise_sim(reps, models, ddf, params, nsim, drop_nonconverged)
  an <- intersect(models, .analytic_models)
  summary$analytical <- NA_real_
  if (length(an) && design$n_placebo > 0) {
    ap <- cqtc_power(design, params, models = an,
                     n_per_cohort = al$n_active + al$n_placebo,
                     n_conc_sim = n_conc_sim, seed = seed)
    summary$analytical <- ap$power[match(summary$model, ap$model)]
  }
  structure(list(summary = summary,
                 replicates = if (keep_replicates) reps else NULL,
                 design = design, params = params, nsim = nsim,
                 n_per_cohort = al$n_active + al$n_placebo),
            class = "cqtc_sim")
}

.summarise_sim <- function(reps, models, ddf, params, nsim, drop_nonconverged = FALSE) {
  U0 <- params$U0
  U1 <- params$U1
  rows <- list()
  for (m in models) {
    r <- if (is.null(reps)) NULL else reps[reps$model == m, , drop = FALSE]
    n_fit <- if (is.null(r)) 0L else nrow(r)
    n_nonconv <- if (n_fit) sum(!r$converged) else 0L
    n_singular <- if (n_fit) sum(r$singular) else 0L
    if (drop_nonconverged && n_fit) r <- r[r$converged, , drop = FALSE]
    n_ok <- if (is.null(r)) 0L else nrow(r)
    for (k in ddf) {
      tag <- .ddf_tag(k)
      if (n_ok > 0) {
        up <- r[[paste0("upper_", tag)]]
        lo <- r[[paste0("lower_", tag)]]
        pw <- mean(up < U0)
        cp <- mean(lo < U1 & up > U1)
        b <- mean(r$estimate - U1)
        rb <- if (U1 != 0) mean((r$estimate - U1) / U1) else NA_real_
      } else {
        pw <- cp <- b <- rb <- NA_real_
      }
      rows[[length(rows) + 1]] <- data.frame(
        model = m, ddf = k, power = pw,
        power_mcse = if (n_ok > 0) sqrt(pw * (1 - pw) / n_ok) else NA_real_,
        coverage = cp, bias = b, rel_bias = rb,
        n_ok = n_ok, n_failed = nsim - n_fit, n_nonconv = n_nonconv,
        n_singular = n_singular, stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Simulation-based sample size search.
.sample_size_sim <- function(design, params, power, models, max_n, n_conc_sim,
                             seed, nsim, ddf, parallel) {
  ddf <- .match_ddf(ddf)
  if (length(ddf) != 1) stop("Give a single `ddf` method for the simulation search.", call. = FALSE)
  b <- design$block
  start_model <- sub("^RS( |$)", "RS0\\1", models)
  start <- cqtc_sample_size(design, params, power = power,
                            models = unique(start_model), method = "analytical",
                            max_n_per_cohort = max_n, n_conc_sim = n_conc_sim,
                            seed = seed)
  rows <- lapply(seq_along(models), function(j) {
    m <- models[j]
    cache <- list()
    emp <- function(n) {
      key <- as.character(n)
      if (is.null(cache[[key]])) {
        s <- cqtc_simulate(design, params, nsim = nsim, models = m, ddf = ddf,
                           n_per_cohort = n, seed = seed, parallel = parallel,
                           keep_replicates = FALSE, n_conc_sim = n_conc_sim)
        cache[[key]] <<- s$summary
      }
      cache[[key]]
    }
    n <- start$n_per_cohort[match(start_model[j], start$model)]
    while (emp(n)$power < power && n + b <= max_n) n <- n + b
    while (n - b >= b && emp(n - b)$power >= power) n <- n - b
    s <- emp(n)
    if (s$power < power) {
      warning(sprintf("%s: target power not reached within max_n_per_cohort = %d.", m, max_n),
              call. = FALSE)
    }
    al <- .allocation(design, n)
    data.frame(model = m, n_active = al$n_active, n_placebo = al$n_placebo,
               n_per_cohort = n, n_total = n * design$n_cohort,
               power = s$power, power_mcse = s$power_mcse, ddf = ddf,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  structure(out, class = c("cqtc_sample_size", "data.frame"),
            design = design, params = params, target = power,
            method = "simulation", nsim = nsim)
}
