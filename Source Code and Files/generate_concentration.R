#' one_compartmental
#'
#' This function generates the concentration data from 1-compartmental model
#'
#' @param dose a scalar of dose level, unit mg
#' @param t a scalar or vector of measurement time
#' @param Ka parameter in one compartmental model
#' @param CL parameter in one compartmental model
#' @param V parameter in one compartmental model
#' @param sd subject variability
#' @param method 1: without log normal error 2: with log normal error
#' @param cv Coefficient of variation of log-normal distribution
#' @return A vector of Concentration
#' @export
one_compartmental <- function(dose, t, Ka = 0.8, CL = 12000, V = 40000, sd = 0.1, method = 1, cv = 0.1){
  eta1 <- rnorm(1, 0, sd)
  eta2 <- rnorm(1, 0, sd)
  eta3 <- rnorm(1, 0, sd)
  Ka <- Ka * exp(eta1)
  CL <- CL * exp(eta2)
  V <- V * exp(eta3)
  Ke <- CL/V
  C <- dose * Ka/(V * (Ka - Ke)) * (exp(-Ke*t) - exp(-Ka*t))
  if(method == 1){
    return(C)
  }
  if(method == 2){
    C <- C * exp(rnorm(length(C), 0, cv))  # how to add log-normal error
  }
  return(C)
}

#' one_compartmental max
#'
#' This function calculates the Cmax of concentration data from 1-compartmental model
#'
#' @param dose a scalar of dose level, unit mg
#' @param t a scalar or vector of measurement time
#' @param Ka parameter in one compartmental model
#' @param CL parameter in one compartmental model
#' @param V parameter in one compartmental model
#' @param sd subject variability
#' @param method 1: without log normal error 2: with log normal error
#' @param cv Coefficient of variation of log-normal distribution
#' @return A vector of Concentration
#' @export
one_compartmental_Cmax <- function(dose = 800, Ka = 0.8, CL = 12000, V = 40000){
  Ke <- CL/V
  Tmax <- (log(Ka) - log(Ke))/(Ka - Ke)
  Cmax <- dose * Ka/(V * (Ka - Ke)) * (exp(-Ke*Tmax) - exp(-Ka*Tmax))
  return(list(Cmax = Cmax, Tmax = Tmax))
}

#' profile_percentage
#'
#' This function generates the concentration data from pre-specified percentage
#'
#' @param Cmax maximum Concentration level, unit mg/mL
#' @param profile pre-specified profile
#' @param cv_C  Coefficient of variation of true concentration
#' @param method 1: without log normal error 2: with log normal error
#' @param cv_error Coefficient of variation of log-normal distribution
#' @return A vector of Concentration
#' @export
profile_percentage <- function(Cmax, profile, cv_C = 0.1, method = 1, cv_error = 0.1){
  C <- Cmax * profile
  for(i in 1:length(C)){
    C[i] <- rnorm(1, C[i], C[i] * cv_C)
  }
  if(method == 1){
    return(C)
  }
  if(method == 2){
    C <- C * exp(rnorm(length(C), 0, cv))
    return(C)
  }

}
