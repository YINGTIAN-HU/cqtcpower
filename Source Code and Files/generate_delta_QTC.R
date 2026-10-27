#' simulate_detla_QTC
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
#' @param C Concentration Profile
#' @return delta QTC and design Matrix
#' @export
simulate_detla_QTC <- function(theta0 = -5.76,
                               theta_C = 1.15/3.33333,
                               theta_Z = 0.5,
                               theta_W = -0.58,
                               theta_T = c(0, 1.93, 2.29, 3.53, -3.15, -0.28, 5.17, 3.98),
                               eta0 = 8,
                               eta1 = 7.5,
                               epsilon = 5,
                               mean_baseline = 408,
                               sd_baseline = 16,
                               rho = 0,
                               C
){

  T <- nrow(C)
  N <- ncol(C)
  delta_QTC <- array(0, dim = c(T, N))
  for(i in 1:N){

    if(sum(C[,i]) == 0){
      # if the subject is in the placebo group
      X <- cbind(1, C[,i], 0, 1, diag(T))
    } else {
      # if the subject is in the treatment group
      X <- cbind(1, C[,i], 1, 1, diag(T))
    }
    QTC_baseline <- rnorm(1, mean_baseline, sd_baseline)
    X[, 4] <- X[, 4] * (QTC_baseline - mean_baseline)


    # generate delta QTC data
    Sigma_MVN <- matrix(c(eta0^2, rho * eta0 * eta1, rho * eta0 * eta1, eta1^2), nrow = 2)
    eta <- mvrnorm(n=1, mu = c(0, 0), Sigma = Sigma_MVN)
    delta_QTC[,i] <- X %*% c(theta0, theta_C, theta_Z, theta_W, theta_T) + eta[1] + eta[2] * X[,2] + rnorm(T, 0, epsilon)

    if(i == 1){
      design <- X
    } else {
      design <- rbind(design, X)
    }
  }
  return(list(delta_QTC = delta_QTC, design = design))
}
