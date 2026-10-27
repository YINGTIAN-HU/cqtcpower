# ---------
# Function used to calcualte sample size
# ---------

calculate_analytical_sample_size2  <- function(N_G = 8, #number of subject in each group
         R = 3, #N_A: N_P
         Cmax = 0.008, #maximum averaged concetration level
         C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
         profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01),
         delta = 5,
         sigma_e = 6,
         sigma0  = 3,
         sigma1 = 600,
         C_supra = 0.008,
         theta0 = 0.5,
         beta = 0.8){
  
  # set up parameter
  Cmax_group <- Cmax * C_alpha
  C <- profile %*% t(Cmax_group) #concenration level in each group
  G <- ncol(C) #number of group/cohort
  T <- nrow(C) #number of observations
  rho <- R/(R+1) #proportion of N_A
  N_A <- N_G * R/(R+1); N_P <- N_G * 1/(R+1)
  
  # get Concentration data
  Concentration <- NULL
  for(i in 1:ncol(C)){
    Concentration <- c(Concentration, rep(C[, i], N_A))
    Concentration <- c(Concentration, rep(rep(0, T), N_P))
  }
  
  # get dataset
  Group <- rep(1:G, each = N_G * T)
  Subject <- rep(1:(N_G*G), each = T)
  Time <- rep(1:T, (N_G*G))
  Baseline <- rep(1, N_G*G*T)
  Treatment <- rep(c(rep(1, N_A * T), rep(0, N_P * T)),  G)
  data <- data.frame(cbind(Concentration, Group, Subject, Time, Baseline, Treatment))
  
  C <- matrix(data$Concentration, nrow = T)
  gi <- calculate_gi(sigma_e = sigma_e, sigma0 = sigma0, sigma1 = sigma1, C = C)
  
  C.subjsum.mean <- mean(apply(C, 2, sum))
  C.subjsum.var <- var(apply(C, 2, sum))
  C.time.var <- apply(C, 1, var)
  
  Cstar.subjsum.mean <- mean(apply(t(t(C)/(1 + gi)), 2, sum))
  Cstar.subjsum.var <- var(apply(t(t(C)/(1 + gi)), 2, sum))
  Cstar.time.var <- apply(t(t(C)/(1 + gi)), 1, var)
  
  (mean_C <- C.subjsum.mean)
  (var_C <- C.subjsum.var)
  (var_CT <- C.time.var)
  
  (mean_C_star <- Cstar.subjsum.mean)
  (var_C_star <- Cstar.subjsum.var)
  (var_CT_star <- Cstar.time.var)
  
  mean_gi2 <- mean(gi^2/(1 + gi)^2)
  
  #G <- length(dose_level)
  
  U0 <- 10
  covariance_intercept_inv <- calculate_variance_random_intercept(sigma_e = sigma_e, sigma0 = sigma0,
                                                                  C = NULL, Z = NULL, W = NULL,
                                                                  var_Z = NULL, var_C = var_C, var_CT = var_CT, mean_C = mean_C,
                                                                  cov_WC = NULL, cov_WZ = NULL, mean_W = NULL, var_W = NULL, N = G * N_G, R = R,
                                                                  simple = TRUE)
  
  ### base Model
  covariance_base <- 1/(covariance_intercept_inv[2,2])
  N_base <- ((qnorm(beta) + qnorm(0.95))/delta)^2 * C_supra^2 * covariance_base * G * N_G
  
  
  ### random Intercept Model
  covariance_intercept <- solve(covariance_intercept_inv)
  N_intercept <- ((qnorm(beta) + qnorm(0.95))/delta)^2 * t(c(1, C_supra)) %*% covariance_intercept %*% c(1, C_supra) * G * N_G
  
  
  covariance_slope_inv <- calculate_variance_random_slope(sigma_e = sigma_e, sigma0 = sigma0, sigma1 = sigma1,
                                                          C = NULL, Z = NULL, W = NULL, mean_gi2 = mean_gi2,
                                                          var_Z = NULL, var_Cstar = var_C_star, var_CstarT = var_CT_star, mean_Cstar = mean_C_star,
                                                          cov_WA = NULL, cov_WCstar = NULL, mean_W = NULL, var_W = NULL, N = G * N_G, R = R,
                                                          simple = TRUE, approx = FALSE)

  
  ### random Slope no treatment Model
  covariance_base_slope <- 1/covariance_slope_inv[2,2]
  N_base_slope <- ((qnorm(beta) + qnorm(0.95))/delta)^2 * C_supra^2 * covariance_base_slope * G * N_G
  
  ### random Slope Model
  covariance_slope <- solve(covariance_slope_inv)
  N_slope <- ((qnorm(beta) + qnorm(0.95))/delta)^2 * t(c(1, C_supra)) %*% covariance_slope %*% c(1, C_supra) * G * N_G
  
  
  N_base <-  ceiling(N_base/(G * (R + 1))) * (R + 1) * G
  N_intercept <- ceiling(N_intercept/(G * (R + 1))) * (R + 1) * G
  N_base_slope <- ceiling(N_base_slope/(G * (R + 1) ))  * (R + 1) * G
  N_slope <- ceiling(N_slope/(G * (R + 1)))  * (R + 1) * G
  
  return(list(N_base = N_base,
              N_intercept = N_intercept,
              N_base_slope = N_base_slope,
              N_slope = N_slope))
}
