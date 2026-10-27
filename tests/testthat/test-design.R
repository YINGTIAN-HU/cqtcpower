test_that("design defaults and blocks", {
  d <- cqtc_design()
  expect_equal(d$n_cohort, 8)
  expect_equal(d$n_time, 7)
  expect_equal(d$block, 4)
  expect_equal(cqtc_design(n_active = 4, n_placebo = 2)$block, 3)
  expect_equal(cqtc_design(n_active = 3, n_placebo = 2)$block, 5)
  expect_equal(cqtc_design(conc_profile())$n_time, 8)
})

test_that("allocation keeps the active:placebo ratio", {
  d <- cqtc_design(n_active = 6, n_placebo = 2)
  expect_equal(cqtcpower:::.allocation(d, 12), list(n_active = 9L, n_placebo = 3L))
  expect_error(cqtcpower:::.allocation(d, 10), "multiple of 4")
})

test_that("input checks", {
  expect_error(cqtc_design(conc_profile(), times = 1:3), "one entry per element")
  expect_error(cqtc_design(n_active = 0), "n_active")
  expect_error(cqtc_params(U1 = 5, rho = 2), "rho")
  expect_error(conc_pk_1cmt(Ka = 0.3, CL = 12000, V = 40000), "differ")
})

test_that("parameters derive theta_C and sigma1", {
  p <- cqtc_params(U1 = 5, theta_Z = 0.5, slope_cv = 0.75)
  expect_equal(p$theta_C, 4.5 / 0.008)
  expect_equal(p$sigma1, 0.75 * 4.5 / 0.008)
  expect_equal(cqtc_params(U1 = 5, sigma1 = 600)$sigma1, 600)
})

test_that("typical concentrations", {
  ct <- conc_typical(cqtc_design())
  expect_equal(dim(ct$conc), c(7, 8))
  expect_gt(ct$true_cmax, max(ct$cmax) * 0.95)
})
