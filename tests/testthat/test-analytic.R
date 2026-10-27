# Exact GLS standard error of U for a given concentration matrix, used as the
# reference the closed-form approximations are checked against.
exact_se_U <- function(C, W, sigma_e, sigma0, sigma1, with_Z, C_supra) {
  T <- nrow(C)
  S <- 0
  for (i in seq_len(ncol(C))) {
    Z <- as.numeric(sum(C[, i]) > 0)
    X <- cbind(if (with_Z) rep(Z, T), C[, i], rep(W[i], T), diag(T))
    V <- sigma_e^2 * diag(T) + sigma0^2 + sigma1^2 * tcrossprod(C[, i])
    S <- S + crossprod(X, solve(V, X))
  }
  Cov <- solve(S)
  if (with_Z) {
    cc <- c(1, C_supra)
    sqrt(drop(t(cc) %*% Cov[1:2, 1:2] %*% cc))
  } else {
    C_supra * sqrt(Cov[1, 1])
  }
}

test_that("closed-form SE matches exact GLS (profile design)", {
  des <- cqtc_design(conc_profile(), 6, 2)
  set.seed(11)
  st <- cqtcpower:::.study_conc(des, 6, 2)
  W <- rnorm(ncol(st$C), 0, 16)
  for (s1 in c(0, 600)) {
    par <- cqtc_params(U1 = 5, sigma1 = s1)
    pw <- cqtc_power(des, par)
    for (m in pw$model) {
      if (s1 == 0 && grepl("RS0", m)) next
      slope <- if (grepl("RS0", m)) s1 else 0
      ex <- exact_se_U(st$C, W, 5, 8, slope, grepl("Z", m), 0.008)
      expect_lt(abs(pw$se_U[pw$model == m] / ex - 1), 0.03)
    }
  }
})

test_that("reference values (paper Figure 2 setting, corrected formula)", {
  des <- cqtc_design(conc_profile(), 6, 2)
  pw <- cqtc_power(des, cqtc_params(U1 = 5, sigma1 = 600))
  expect_equal(pw$se_U, c(1.813641611, 2.795178543, 2.180697134, 3.044301620),
               tolerance = 1e-7)
  # closed-form cohort sizes for ddQTc = 0..5 (delta = 10..5), as produced by
  # calculate_analytical_sample_size2() with C/(1 + g_i)
  ref <- rbind(c(4, 4, 4, 8), c(4, 8, 4, 8), c(4, 8, 4, 8),
               c(4, 8, 8, 12), c(8, 12, 8, 16), c(8, 16, 12, 20))
  for (k in 0:5) {
    ss <- cqtc_sample_size(des, cqtc_params(U1 = k, sigma1 = 600))
    expect_equal(ss$n_formula, ref[k + 1, ])
  }
})

test_that("sample size reaches target and is minimal", {
  des <- cqtc_design(conc_profile(), 6, 2)
  par <- cqtc_params(U1 = 4, sigma1 = 600)
  ss <- cqtc_sample_size(des, par, power = 0.9)
  expect_true(all(ss$power >= 0.9))
  for (i in seq_len(nrow(ss))) {
    if (ss$n_per_cohort[i] > des$block) {
      smaller <- cqtc_power(des, par, models = ss$model[i],
                            n_per_cohort = ss$n_per_cohort[i] - des$block)
      expect_lt(smaller$power, 0.9)
    }
  }
  expect_error(cqtc_sample_size(des, cqtc_params(U1 = 10)), "below the threshold")
})

test_that("power curve is monotone and models are ordered", {
  pc <- cqtc_power_curve(cqtc_design(conc_profile()), cqtc_params(U1 = 5, sigma1 = 600))
  for (m in unique(pc$model)) expect_true(all(diff(pc$power[pc$model == m]) >= 0))
  at8 <- pc[pc$n_per_cohort == 8, ]
  expect_true(at8$power[at8$model == "RI"] >= at8$power[at8$model == "RS0"])
  expect_true(at8$power[at8$model == "RI + Z"] >= at8$power[at8$model == "RS0 + Z"])
  expect_silent(plot(pc))
})

test_that("type I error is alpha at U1 = U0", {
  pw <- cqtc_power(cqtc_design(conc_profile()), cqtc_params(U1 = 10, sigma1 = 600))
  expect_equal(pw$power, rep(0.05, 4), tolerance = 1e-12)
})

test_that("random PK concentrations: seed gives reproducible moments and restores RNG", {
  des <- cqtc_design(conc_pk_1cmt(), 6, 2)
  par <- cqtc_params(U1 = 5, theta_Z = 0.5, slope_cv = 0.75)
  set.seed(99)
  before <- .Random.seed
  a <- cqtc_power(des, par, n_conc_sim = 50, seed = 3)
  expect_identical(.Random.seed, before)
  b <- cqtc_power(des, par, n_conc_sim = 50, seed = 3)
  expect_equal(a$power, b$power)
})

test_that("mean_gi2 averages over replicates (bug fix)", {
  des <- cqtc_design(conc_pk_1cmt(), 6, 2)
  par <- cqtc_params(U1 = 5, theta_Z = 0.5, slope_cv = 0.75)
  m <- conc_moments(des, par, n_conc_sim = 20, seed = 1)
  one <- cqtcpower:::.with_seed(1, {
    vals <- vapply(1:20, function(i) {
      st <- cqtcpower:::.study_conc(des, 6, 2)
      cqtcpower:::.conc_moments_one(st$C, 5, 8, par$sigma1)$mean_gi2
    }, numeric(1))
  })
  expect_equal(m$mean_gi2, mean(vals))
})

test_that("RS models are refused analytically", {
  expect_error(cqtc_power(cqtc_design(), cqtc_params(U1 = 5), models = "RS + Z"),
               "cqtc_simulate")
  expect_equal(cqtcpower:::.match_models(c("ri+z", "RS0 + Z")), c("RI + Z", "RS0 + Z"))
})
