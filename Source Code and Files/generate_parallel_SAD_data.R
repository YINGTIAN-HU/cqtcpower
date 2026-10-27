#' parallel_SAD_sim
#'
#' This function generates the detla QTC directly based on the profile data C.
#'
#' @param theta_0 population mean intercept in the absence of the treatment effect
#' @param theta_C population mean slope of the assumed linear association between concentration and delta QTC
#' @param theta_Z fixed effect associated with treatment
#' @param theta_W fixed effect associated with baseline
#' @param theta_T fixed effect associated with time
#' @param eta_0 random intercept
#' @param eta_1 random slope
#' @param epsilon random error
#' @param mean_baseline mean of baseline Concentration
#' @param sd_baseline sd of baseline Concentration
#' @param N_A Number of treatment subject in each cohort
#' @param N_P Number of placebo subject in each cohort
#' @param t ECG measurment time points
#' @param dose_level dose level of each cohort
#' @param C_supra supratherapeutic Concentration
#' @param fit_model True: fit LMM; false: just generate data not fit model
#' @param model intercept: random intercept model; slope: random slope model
#' @return delta QTC and design Matrix
#' @export
# t <- c(1, 2, 4, 8, 12, 16, 24)
# T <- length(t)
# N_A <- 12
# N_P <- 4
# dose_level <- c(1, 3, 10, 30, 100, 200, 400, 800)
# C_supra <- 0.008

parallel_SAD_sim <- function(theta0 = -5.76, #
                             theta_C = 1.15/3.33333,
                             theta_Z = 0.5,
                             theta_W = -0.58,
                             theta_T = c(1.93, 2.29, 3.53, -3.15, -0.28, 5.17, 3.98),
                             eta0 = 8,
                             eta1 = 7.5,
                             epsilon = 5,
                             mean_baseline = 408,
                             sd_baseline = 16,
                             N_A = 6,
                             N_P = 2,
                             t = c(1, 2, 4, 8, 12, 16, 24),
                             dose_level = c(1, 3, 10, 30, 100, 200, 400, 800),
                             C_supra = 0.008,
                             rho = 0,
                             fit_model = TRUE,
                             model = "slope",
                             Ka = 0.8, CL = 12000, V = 40000, sd = 0.1, method = 1, cv = 0.1){
  delta_QTC <- list()
  design <- list()
  cohort_num <- length(dose_level)
  T <- length(t)
  for(cohort in 1:cohort_num){

    # generate treatment concentration data
    C <- array(0, dim = c(T, N_A))
    for(i in 1:N_A){
      C[, i] <- one_compartmental(dose_level[cohort], t, Ka = Ka, CL = CL, V = V, sd = sd, method = method, cv = cv)
    }
    # generate placebo concentration data
    C <- cbind(C, array(0, dim = c(T, N_P)))

    # generate Delta QTC data
    simdata <- simulate_detla_QTC(theta0 = theta0,
                                  theta_Z = theta_Z,
                                  theta_C = theta_C,
                                  theta_W = theta_W,
                                  theta_T = theta_T,
                                  eta0 = eta0,
                                  eta1 = eta1,
                                  epsilon = epsilon,
                                  mean_baseline = mean_baseline,
                                  sd_baseline = sd_baseline,
                                  C = C, rho = rho)

    # get each cohort data
    delta_QTC[[cohort]] <- c(simdata$delta_QTC)
    design[[cohort]] <- cbind(simdata$design, dose_level[cohort])
  }
  design <- do.call(rbind, design)
  delta_QTC <- unlist(delta_QTC)

  data <- data.frame(cbind(delta_QTC, design))
  colnames(data) <- c("delta_QTC", "Intercept", "Concentration", "Treatment", "Baseline",  paste("t", 1:T, sep = ""), "Dose")
  data$Subject <- rep(1:(cohort_num *(N_A + N_P)), each = T)
  data <- data[, -2] # remove intercept column

  if(fit_model == TRUE){
    if(model == "RI"){
      # formula <- as.formula(paste("delta_QTC ~ -1 + Concentration + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "), sep = ""))
      # fit <- lme(formula, data = data, random = ~ 1|Subject, control = lmeControl(opt = "optim"))
      formula <- as.formula(paste("delta_QTC ~ -1 + Concentration + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "),  " + (1|Subject)", sep = ""))
      fit <- lmer(formula, data = data, control = lmerControl(optimizer ="bobyqa"))
    }
    if(model == "RI + Z"){
      #formula <- as.formula(paste("delta_QTC ~ -1 + Treatment + Concentration + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "), sep = ""))
      #fit <- lme(formula, data = data, random = ~ 1|Subject, control = lmeControl(opt = "optim"))
      formula <- as.formula(paste("delta_QTC ~ -1 + Treatment + Concentration + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "),  " + (1|Subject)", sep = ""))
      fit <- lmer(formula, data = data, control = lmerControl(optimizer ="bobyqa"))
    }
    if(model == "RS0"){
      formula <- as.formula(paste("delta_QTC ~ -1 + Concentration + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "),  "+ (Concentration||Subject)", sep = ""))
      fit <- lmer(formula, data = data, control = lmerControl(optimizer ="bobyqa"))
    }
    if(model == "RS0 + Z"){
      formula <- as.formula(paste("delta_QTC ~ -1 + Treatment + Concentration + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "),  " + (Concentration||Subject)", sep = ""))
      fit <- lmer(formula, data = data, control = lmerControl(optimizer ="bobyqa"))
    }
    if(model == "RS"){
      formula <- as.formula(paste("delta_QTC ~ -1 + Concentration + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "),  " + (Concentration|Subject)", sep = ""))
      fit <- lmer(formula, data = data, control = lmerControl(optimizer ="bobyqa"))
    }
    if(model == "RS + Z"){
      #formula <- as.formula(paste("delta_QTC ~ -1 + Treatment + Concentration + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "), sep = ""))
      #fit <- lme(formula, data = data, random = ~ 1 + Concentration|Subject, control = lmeControl(opt = "optim"))
      formula <- as.formula(paste("delta_QTC ~ -1 + Treatment + Concentration + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "),  " + (Concentration|Subject)", sep = ""))
      fit <- lmer(formula, data = data, control = lmerControl(optimizer ="bobyqa"))
    }


    if(model != "RI" & model != "RS0" & model != "RS"){
      Sigma <- vcov(fit)[c(1,2), c(1,2)]
      #if(model == "slope0"){
      coef <- summary(fit)$coefficients
      #} else {
      #  coef <- summary(fit)$coefficients$fixed
      #}
      est <- coef[1] + C_supra * coef[2]
      std <- sqrt(Sigma[1, 1] + C_supra^2 * Sigma[2,2] + 2*C_supra * Sigma[1, 2])
      CI <- c(est - qnorm(0.95) *std, est + qnorm(0.95) *std)

      KRDF <- unlist(contest(fit, c(1, C_supra, 0, rep(0, T)), ddf = "Kenward-Roger")[4])
      SWDF <- unlist(contest(fit, c(1, C_supra, 0, rep(0, T)))[4])
      CIKR <- c(est - qt(0.95, KRDF) *std, est +  qt(0.95, KRDF) *std)
      CISW <- c(est - qt(0.95, SWDF) *std, est +  qt(0.95, SWDF) *std)

    } else {
      Sigma <- vcov(fit)[1,1]
      #if(model == "base_slope"){
      coef <- summary(fit)$coefficients
      #} else {
      #  coef <- summary(fit)$coefficients$fixed
      #}
      est <- C_supra * coef[1]
      std <- C_supra * sqrt(Sigma)
      CI <- c(est - qnorm(0.95) *std, est + qnorm(0.95) *std)

      KRDF <- unlist(contest(fit, c(C_supra, 0, rep(0, T)), ddf = "Kenward-Roger")[4])
      SWDF <- unlist(contest(fit, c(C_supra, 0, rep(0, T)))[4])
      CIKR <- c(est - qt(0.95, KRDF) *std, est +  qt(0.95, KRDF) *std)
      CISW <- c(est - qt(0.95, SWDF) *std, est +  qt(0.95, SWDF) *std)
    }
    res <- data.frame("est" = est, "std" = std,
                      "lower" = CI[1], "upper" = CI[2],
                      "lowerKR" = CIKR[1], "upperKR" = CIKR[2],
                      "lowerSW" = CISW[1], "upperSW" = CISW[2])
    return(list(data = data, res = res, fit = fit))
  } else {
    return(list(data = data))
  }
}
