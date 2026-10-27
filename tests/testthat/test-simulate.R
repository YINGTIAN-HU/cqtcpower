test_that("simulated data have the expected structure", {
  des <- cqtc_design()
  d <- simulate_sad_data(des, cqtc_params(U1 = 5, theta_Z = 0.5), seed = 1)
  expect_equal(nrow(d), 8 * 8 * 7)
  expect_equal(nlevels(d$Subject), 64)
  expect_true(all(d$Concentration[d$Treatment == 0] == 0))
  expect_true(all(d$Concentration[d$Treatment == 1] > 0))
  expect_identical(d, simulate_sad_data(des, cqtc_params(U1 = 5, theta_Z = 0.5), seed = 1))
  d12 <- simulate_sad_data(des, cqtc_params(U1 = 5), n_per_cohort = 12, seed = 2)
  expect_equal(nlevels(d12$Subject), 96)
})

test_that("data follow the generating model (noise-free check)", {
  des <- cqtc_design(conc_pk_1cmt(iiv_sd = 0))
  par <- cqtc_params(U1 = 6, theta_Z = 1, sigma_e = 1e-8, sigma0 = 0,
                     baseline_sd = 0, theta0 = 0, theta_T = rep(0, 7))
  d <- simulate_sad_data(des, par, seed = 1)
  expect_equal(d$delta_QTC, par$theta_Z * d$Treatment + par$theta_C * d$Concentration,
               tolerance = 1e-6)
})

test_that("model fits return estimates and bounds", {
  d <- simulate_sad_data(cqtc_design(), cqtc_params(U1 = 5, theta_Z = 0.5), seed = 4)
  f <- fit_cqtc_model(d, "RI + Z", C_supra = 0.008, ddf = c("Wald", "Satterthwaite"))
  expect_named(f, c("model", "estimate", "se", "converged", "singular", "optimizer",
                    "df_Wald", "lower_Wald", "upper_Wald",
                    "df_SW", "lower_SW", "upper_SW"))
  expect_true(f$upper_SW > f$upper_Wald)
  expect_equal(f$upper_Wald - f$estimate, qnorm(0.95) * f$se)
})

test_that("C_rel scaling is an exact reparameterisation", {
  d <- simulate_sad_data(cqtc_design(), cqtc_params(U1 = 5, theta_Z = 0.5), seed = 4)
  f <- fit_cqtc_model(d, "RI + Z", C_supra = 0.008, ddf = "Wald")
  raw <- suppressWarnings(lme4::lmer(
    delta_QTC ~ 0 + Time + Baseline + Treatment + Concentration + (1 | Subject),
    data = d, control = lme4::lmerControl(optimizer = "bobyqa")))
  L <- c(Treatment = 1, Concentration = 0.008)
  b <- lme4::fixef(raw)[names(L)]
  V <- as.matrix(stats::vcov(raw))[names(L), names(L)]
  expect_equal(f$estimate, sum(L * b), tolerance = 1e-6)
  expect_equal(f$se, sqrt(drop(t(L) %*% V %*% L)), tolerance = 1e-5)
})

test_that("random-slope fits converge on the paper's design", {
  des <- cqtc_design()
  par <- cqtc_params(U1 = 5, theta_Z = 0.5, slope_cv = 0.75)
  fits <- do.call(rbind, lapply(1:8, function(i) {
    fit_cqtc_model(simulate_sad_data(des, par, seed = i), "RS + Z", 0.008, ddf = "Wald")
  }))
  expect_gte(mean(fits$converged), 0.875)
})

test_that("cqtc_simulate summarises and is reproducible", {
  des <- cqtc_design()
  par <- cqtc_params(U1 = 5, theta_Z = 0.5)
  s1 <- cqtc_simulate(des, par, nsim = 4, models = c("RI", "RI + Z"), ddf = "Wald",
                      n_conc_sim = 10)
  s2 <- cqtc_simulate(des, par, nsim = 4, models = c("RI", "RI + Z"), ddf = "Wald",
                      n_conc_sim = 10)
  expect_s3_class(s1, "cqtc_sim")
  expect_equal(s1$summary, s2$summary)
  expect_equal(nrow(s1$summary), 2)
  expect_true(all(s1$summary$n_ok == 4))
  expect_false(anyNA(s1$summary$analytical))
  expect_output(print(s1), "Simulated C-QTc")
})

test_that("ddf names are normalised", {
  expect_equal(cqtcpower:::.match_ddf(c("KR", "sw", "wald")),
               c("Kenward-Roger", "Satterthwaite", "Wald"))
  expect_error(cqtcpower:::.match_ddf("foo"), "ddf")
})
