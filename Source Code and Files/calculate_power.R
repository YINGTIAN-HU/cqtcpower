# ------------------
# Functions used to calculate analytic power and empirical power
# ------------------

SummaryConcentration <- function(U1 = 5,
                                 CV = 0.75, # CV of random slope
                                 C_supra = 0.008,
                                 theta_Z_true = 0.5,
                                 theta_T_true = c(1.93, 2.29, 3.53, -3.15, -0.28, 5.17, 3.98),
                                 N_A = 6,
                                 N_P = 2,
                                 t = c(1, 2, 4, 8, 12, 16, 24),
                                 dose_level = c(1, 3, 10, 30, 100, 200, 400, 800),
                                 eta0 = 8,
                                 epsilon = 5, rho =0,  N_sim = 1000,
                                 Ka = 0.8, CL = 12000, V = 40000, sd = 0.1, method = 1, cv = 0.1){

  theta_C_true = (U1 - theta_Z_true)/C_supra # treatment specific effect
  eta1 = CV * theta_C_true

  mean_C_all <- rep(0, N_sim)
  var_C_all <- rep(0, N_sim)
  var_CT_all <- array(0, dim = c(N_sim, length(t)))
  mean_C_star_all <- rep(0, N_sim)
  var_C_star_all <- rep(0, N_sim)
  var_CT_star_all <- array(0, dim = c(N_sim, length(t)))
  mean_gi2_all <- rep(0, N_sim)

  for(n_sim in 1:N_sim){
    set.seed(n_sim)
    res <- parallel_SAD_sim(theta0 = -5.76, #
                            theta_C = theta_C_true,
                            theta_Z = theta_Z_true,
                            theta_W = -0.58,
                            theta_T = theta_T_true,
                            eta0 = eta0,
                            eta1 = eta1,
                            epsilon = epsilon,
                            mean_baseline = 408,
                            sd_baseline = 16,
                            N_A = N_A,
                            N_P = N_P,
                            t = t,
                            dose_level = dose_level,
                            C_supra = C_supra,
                            fit_model = FALSE,
                            model = "intercept", rho = rho,
                            Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv)

    C <- matrix(res$data$Concentration, nrow = length(t))
    gi <- calculate_gi(sigma_e = epsilon, sigma0 = eta0, sigma1 = CV* theta_C_true, C = C)

    C.subjsum.mean <- mean(apply(C, 2, sum))
    C.subjsum.var <- var(apply(C, 2, sum))
    C.time.var <- apply(C, 1, var)

    Cstar.subjsum.mean <- mean(apply(t(t(C)/(1 + gi)), 2, sum))
    Cstar.subjsum.var <- var(apply(t(t(C)/(1 + gi)), 2, sum))
    Cstar.time.var <- apply(t(t(C)/(1 + gi)), 1, var)

    
    mean_C_all[n_sim] <- C.subjsum.mean
    var_C_all[n_sim] <- C.subjsum.var
    var_CT_all[n_sim, ] <- C.time.var

    mean_C_star_all[n_sim] <- Cstar.subjsum.mean
    var_C_star_all[n_sim] <- Cstar.subjsum.var
    var_CT_star_all[n_sim, ] <- Cstar.time.var
    
    mean_gi2_all[n_sim] <- mean(gi^2/(1 + gi)^2)
  }
  mean_C <- mean(mean_C_all)
  var_C <- mean(var_C_all)
  var_CT <- apply(var_CT_all, 2, mean)

  mean_C_star <- mean(mean_C_star_all)
  var_C_star <- mean(var_C_star_all)
  var_CT_star <- apply(var_CT_star_all, 2, mean)
  
  mean_gi2 <- mean(mean_gi2_all)
  
  return(list(mean_C = mean_C,
              var_C = var_C,
              var_CT = var_CT,
              mean_C_star = mean_C_star,
              var_C_star = var_C_star,
              var_CT_star = var_CT_star,
              mean_gi2 = mean_gi2))
}


calculate_analytical_power <- function(U1 = 5,
                                       CV = 0.75, # CV of random slope
                                       C_supra = 0.008,
                                       theta_Z_true = 0.5,
                                       theta_T_true = c(1.93, 2.29, 3.53, -3.15, -0.28, 5.17, 3.98),
                                       N_A = 6,
                                       N_P = 2,
                                       t = c(1, 2, 4, 8, 12, 16, 24),
                                       dose_level = c(1, 3, 10, 30, 100, 200, 400, 800),
                                       eta0 = 8,
                                       epsilon = 5, rho = 0,
                                       Ka = 0.8, CL = 12000, V = 40000, sd = 0.1, method = 1, cv = 0.1){

  theta_C_true = (U1 - theta_Z_true)/C_supra # treatment specific effect
  eta1 = CV * theta_C_true
  concentration_stats <- SummaryConcentration(U1 = U1,
                                              CV = CV, # CV of random slope
                                              C_supra = C_supra,
                                              theta_Z_true = theta_Z_true,
                                              theta_T_true = theta_T_true,
                                              N_A = N_A ,
                                              N_P = N_P ,
                                              t = t,
                                              dose_level = dose_level,
                                              eta0 = eta0,
                                              epsilon = epsilon, rho = rho, N_sim = 1000,
                                              Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv)
  mean_C = concentration_stats$mean_C
  var_C = concentration_stats$var_C
  var_CT = concentration_stats$var_CT
  mean_C_star = concentration_stats$mean_C_star
  var_C_star = concentration_stats$var_C_star
  var_CT_star = concentration_stats$var_CT_star
  mean_gi2 = concentration_stats$mean_gi2

  U0 <- 10
  covariance_intercept_inv <- calculate_variance_random_intercept(sigma_e = epsilon, sigma0 = eta0,
                                                                  C = NULL, Z = NULL, W = NULL,
                                                                  var_Z = NULL, var_C = var_C, var_CT = var_CT, mean_C = mean_C,
                                                                  cov_WC = NULL, cov_WZ = NULL, mean_W = NULL, var_W = NULL, N = length(dose_level) * (N_A + N_P), R = N_A/N_P,
                                                                  simple = TRUE)

  ### base Model
  covariance_base <- 1/(covariance_intercept_inv[2,2])
  variance_base <- C_supra^2 * covariance_base
  power_analytic_base <- pnorm(-qnorm(0.95) + (U0 - U1)/sqrt(variance_base))


  ### random Intercept Model
  covariance_intercept <- solve(covariance_intercept_inv)
  variance_intercept <- t(c(1, C_supra)) %*% covariance_intercept %*% c(1, C_supra)
  power_analytic_intercept <- pnorm(-qnorm(0.95) + (U0 - U1)/sqrt(variance_intercept))



  covariance_slope_inv <- calculate_variance_random_slope(sigma_e = epsilon, sigma0 = eta0, sigma1 = CV * theta_C_true,
                                                          C = NULL, Z = NULL, W = NULL, mean_gi2 = mean_gi2, 
                                                          var_Z = NULL, var_Cstar = var_C_star, var_CstarT = var_CT_star, mean_Cstar = mean_C_star,
                                                          cov_WA = NULL, cov_WCstar = NULL, mean_W = NULL, var_W = NULL, N = length(dose_level) * (N_A + N_P), R = N_A/N_P,
                                                          simple = TRUE, approx = FALSE)

  ### random Slope no treatment Model
  covariance_base_slope <- 1/covariance_slope_inv[2,2]
  variance_base_slope <- C_supra^2 * covariance_base_slope
  power_analytic_base_slope <- pnorm(-qnorm(0.95) + (U0 - U1)/sqrt(variance_base_slope))

  ### random Slope Model
  covariance_slope <- solve(covariance_slope_inv)
  variance_slope <- t(c(1, C_supra)) %*% covariance_slope %*% c(1, C_supra)
  power_analytic_slope <- pnorm(-qnorm(0.95) + (U0 - U1)/sqrt(variance_slope))

  return(list(base = power_analytic_base,
              intercept = power_analytic_intercept,
              base_slope=power_analytic_base_slope,
              slope=power_analytic_slope))
}



calculate_empirical_power <- function(U1 = 5,
                                      CV = 0.75, # CV of random slope
                                      C_supra = 0.008,
                                      theta_Z_true = 0.5,
                                      theta_T_true = c(1.93, 2.29, 3.53, -3.15, -0.28, 5.17, 3.98),
                                      N_A = 6,
                                      N_P = 2,
                                      t = c(1, 2, 4, 8, 12, 16, 24),
                                      dose_level = c(1, 3, 10, 30, 100, 200, 400, 800),
                                      eta0 = 8,
                                      epsilon = 5, N_sim = 1000, rho = 0, model = "base",
                                      Ka = 0.8, CL = 12000, V = 40000, sd = 0.1, method = 1, cv = 0.1){

  theta_C_true = (U1 - theta_Z_true)/C_supra # treatment specific effect
  eta1 = CV * theta_C_true

  res.comb <- NULL
  coef.comb <- NULL
  var.comb <- NULL
  for(i in 1:N_sim){
    # if(i %% 10 == 0){
    #  print(i)
    # }
    set.seed(i)
    res <- tryCatch(parallel_SAD_sim(theta0 = -5.76, #
                                     theta_C = theta_C_true,
                                     theta_Z = theta_Z_true,
                                     theta_W = -0.58,
                                     theta_T = theta_T_true,
                                     eta0 = eta0,
                                     eta1 = eta1,
                                     epsilon = epsilon,
                                     mean_baseline = 408,
                                     sd_baseline = 16,
                                     N_A = N_A,
                                     N_P = N_P,
                                     t = t,
                                     dose_level = dose_level,
                                     C_supra = C_supra,
                                     rho = rho,
                                     fit_model = TRUE,
                                     model = model,
                                     Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv), error = function(c) return(NA))
    if(length(res) > 1){
      #if(model != "slope0" & model != "base_slope"){
      #  coef <- summary(res$fit)$coefficients$fixed
      #} else {
      coef <- summary(res$fit)$coefficients
      #}
      coef.comb <- rbind(coef.comb, coef)
      res.comb <- rbind(res.comb, res$res)
      #var.comb  <- rbind(var.comb,  as.numeric(VarCorr(res$fit)[,2]))
    }
  }
  power <- mean(res.comb$upper <10, na.rm =TRUE)  # power
  power.sd <- sd(res.comb$upper <10, na.rm =TRUE)  # power sd
  CP <- mean(res.comb$upper > U1 & res.comb$lower < U1, na.rm = TRUE) # CP
  CP.sd <- sd(res.comb$upper > U1 & res.comb$lower < U1, na.rm = TRUE) # CP

  power.KR <- mean(res.comb$upperKR <10, na.rm =TRUE)  # power
  power.sd.KR <- sd(res.comb$upperKR <10, na.rm =TRUE)  # power sd
  CP.KR <- mean(res.comb$upperKR > U1 & res.comb$lowerKR < U1, na.rm = TRUE) # CP
  CP.sd.KR <- sd(res.comb$upperKR > U1 & res.comb$lowerKR < U1, na.rm = TRUE) # CP

  power.SW <- mean(res.comb$upperSW <10, na.rm =TRUE)  # power
  power.sd.SW <- sd(res.comb$upperSW <10, na.rm =TRUE)  # power sd
  CP.SW <- mean(res.comb$upperSW > U1 & res.comb$lowerSW < U1, na.rm = TRUE) # CP
  CP.sd.SW <- sd(res.comb$upperSW > U1 & res.comb$lowerSW < U1, na.rm = TRUE) # CP

  bias <- mean((res.comb$est - U1)/U1, na.rm = TRUE) # relative bias
  bias.sd <- sd((res.comb$est - U1)/U1, na.rm = TRUE) # relative bias

  return(list(bias = bias, CP = CP , power = power,
              bias.sd = bias.sd, CP.sd = CP.sd , power.sd = power.sd,
              power.KR = power.KR, power.sd.KR = power.sd.KR,
              CP.KR = CP.KR, CP.sd.KR = CP.sd.KR,
              power.SW = power.SW, power.sd.SW = power.sd.SW,
              CP.SW = CP.SW, CP.sd.SW = CP.sd.SW))
}

calculate_power <- function(U1 = 5,
                            CV = 0, # CV of random slope
                            theta_Z_true = 0,
                            C_supra = 0.008,
                            N_A = 6,
                            N_P = 2,
                            t = c(1, 2, 4, 8, 12, 16, 24),
                            dose_level = c(1, 3, 10, 30, 100, 200, 400, 800),
                            eta0 = 8,
                            epsilon = 5,
                            rho = 0,
                            theta_T_true = c(1.93, 2.29, 3.53, -3.15, -0.28, 5.17, 3.98),
                            N_sim = 1000,
                            Ka = 0.8, CL = 12000, V = 40000, sd = 0.1, method = 1, cv = 0.1){


  analytical_power <-  calculate_analytical_power(U1 = U1,
                                                  CV = CV, # CV of random slope
                                                  C_supra = C_supra,
                                                  theta_Z_true = theta_Z_true,
                                                  theta_T_true = theta_T_true,
                                                  N_A = N_A,
                                                  N_P = N_P,
                                                  t = t,
                                                  dose_level = dose_level,
                                                  eta0 = eta0,
                                                  epsilon = epsilon,
                                                  rho = rho,
						  Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv)


  empirical_power_RI  <- calculate_empirical_power(U1 = U1,
                                                   CV = CV, # CV of random slope
                                                   C_supra = C_supra,
                                                   theta_Z_true = theta_Z_true,
                                                   theta_T_true = theta_T_true,
                                                   N_A = N_A,
                                                   N_P = N_P,
                                                   t = t,
                                                   dose_level = dose_level,
                                                   eta0 = eta0,
                                                   epsilon = epsilon,
                                                   rho = rho,
                                                   N_sim = N_sim, model = "RI",
                                                   Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv)

  empirical_power_RI_Z  <- calculate_empirical_power(U1 = U1,
                                                     CV = CV, # CV of random slope
                                                     C_supra = C_supra,
                                                     theta_Z_true = theta_Z_true,
                                                     theta_T_true = theta_T_true,
                                                     N_A = N_A,
                                                     N_P = N_P,
                                                     t = t,
                                                     dose_level = dose_level,
                                                     eta0 = eta0,
                                                     epsilon = epsilon,
                                                     rho = rho,
                                                     N_sim = N_sim, model = "RI + Z",
                                                     Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv)

  empirical_power_RS0  <- calculate_empirical_power(U1 = U1,
                                                    CV = CV, # CV of random slope
                                                    C_supra = C_supra,
                                                    theta_Z_true = theta_Z_true,
                                                    theta_T_true = theta_T_true,
                                                    N_A = N_A,
                                                    N_P = N_P,
                                                    t = t,
                                                    dose_level = dose_level,
                                                    eta0 = eta0,
                                                    epsilon = epsilon,
                                                    rho = rho,
                                                    N_sim = N_sim, model = "RS0",
                                                    Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv)

  empirical_power_RS0_Z  <- calculate_empirical_power(U1 = U1,
                                                      CV = CV, # CV of random slope
                                                      C_supra = C_supra,
                                                      theta_Z_true = theta_Z_true,
                                                      theta_T_true = theta_T_true,
                                                      N_A = N_A,
                                                      N_P = N_P,
                                                      t = t,
                                                      dose_level = dose_level,
                                                      eta0 = eta0,
                                                      epsilon = epsilon,
                                                      rho = rho,
                                                      N_sim = N_sim, model = "RS0 + Z",
                                                      Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv)

  empirical_power_RS  <- calculate_empirical_power(U1 = U1,
                                                   CV = CV, # CV of random slope
                                                   C_supra = C_supra,
                                                   theta_Z_true = theta_Z_true,
                                                   theta_T_true = theta_T_true,
                                                   N_A = N_A,
                                                   N_P = N_P,
                                                   t = t,
                                                   dose_level = dose_level,
                                                   eta0 = eta0,
                                                   epsilon = epsilon,
                                                   rho = rho,
                                                   N_sim = N_sim, model = "RS",
                                                   Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv)

  empirical_power_RS_Z  <- calculate_empirical_power(U1 = U1,
                                                     CV = CV, # CV of random slope
                                                     C_supra = C_supra,
                                                     theta_Z_true = theta_Z_true,
                                                     theta_T_true = theta_T_true,
                                                     N_A = N_A,
                                                     N_P = N_P,
                                                     t = t,
                                                     dose_level = dose_level,
                                                     eta0 = eta0,
                                                     epsilon = epsilon,
                                                     rho = rho,
                                                     N_sim = N_sim, model = "RS + Z",
                                                     Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv)

  power <- data.frame(analytical = c(analytical_power$base,
                                     analytical_power$intercept,
                                     analytical_power$base_slope,
                                     analytical_power$slope, NA, NA),
                      empirical = c(empirical_power_RI$power,
                                    empirical_power_RI_Z$power,
                                    empirical_power_RS0$power,
                                    empirical_power_RS0_Z$power,
                                    empirical_power_RS$power,
                                    empirical_power_RS_Z$power),
                      empiricalKR = c(empirical_power_RI$power.KR,
                                    empirical_power_RI_Z$power.KR,
                                    empirical_power_RS0$power.KR,
                                    empirical_power_RS0_Z$power.KR,
                                    empirical_power_RS$power.KR,
                                    empirical_power_RS_Z$power.KR),
                      empiricalSW = c(empirical_power_RI$power.SW,
                                    empirical_power_RI_Z$power.SW,
                                    empirical_power_RS0$power.SW,
                                    empirical_power_RS0_Z$power.SW,
                                    empirical_power_RS$power.SW,
                                    empirical_power_RS_Z$power.SW))

  CP <- data.frame(CP = c(empirical_power_RI$CP,
                          empirical_power_RI_Z$CP,
                          empirical_power_RS0$CP,
                          empirical_power_RS0_Z$CP,
                          empirical_power_RS$CP,
                          empirical_power_RS_Z$CP),
                   CPKR = c(empirical_power_RI$CP.KR,
                          empirical_power_RI_Z$CP.KR,
                          empirical_power_RS0$CP.KR,
                          empirical_power_RS0_Z$CP.KR,
                          empirical_power_RS$CP.KR,
                          empirical_power_RS_Z$CP.KR),
                   CPSW = c(empirical_power_RI$CP.SW,
                            empirical_power_RI_Z$CP.SW,
                            empirical_power_RS0$CP.SW,
                            empirical_power_RS0_Z$CP.SW,
                            empirical_power_RS$CP.SW,
                            empirical_power_RS_Z$CP.SW))

  bias <- c(empirical_power_RI$bias,
            empirical_power_RI_Z$bias,
            empirical_power_RS0$bias,
            empirical_power_RS0_Z$bias,
            empirical_power_RS$bias,
            empirical_power_RS_Z$bias)

  return(list(bias = bias, CP = CP, power = power))
}
