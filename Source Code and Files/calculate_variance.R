#' calculate_variance_random_intercept
#'
#' This function is used to calculate the variance estimator of (theta_Z, theta_C) in the random intercept model.
#' There are several senario we want to calculate variance:
#' 1. Known C, Z, W: do not need to specific later parameters. Can choose simple to be TRUE or FALSE
#' 2. Known C, Z: do not need to specify later parameters. Can choose simple to be TRUE
#' 3. Without known C, Z, W: need to specify later variance and covariance parameters and also sample size N. Can choose simple to be TRUE without known information about W
#'
#'
#' @param sigma_e random error of residual
#' @param sigma_0 random intercept variability
#' @param C fixed effect associated with baseline
#' @param Z random intercept
#' @param W random slope
#' @param var_Z variance of treatment Z
#' @param var_C variance of (sum of Concentration level of each subject)
#' @param var_CT vector of variance of Concentration level at each time point
#' @param mean_C mean of (sum of Concentration level of each subject)
#' @param cov_WC covariance of (sum of Concentration level of each subject) and Baseline
#' @param cov_WZ covariance of treatment and Baseline
#' @param mean_W mean of baseline
#' @param var_W variance of baseline
#' @param N Number of subjects
#' @param R Ratio of treatment subjects v.s placebo subjects
#' @param simple TRUE: ignore correlation between treatment and baseline and between sum of Concentration level and Baseline
#' @return variance estimator of (theta_Z, theta_C)
#' @examples
#' var <- calculate_variance_random_intercept(sigma_e, sigma0, C, Z, W, simple = TRUE)
#' var
#' @export
calculate_variance_random_intercept <- function(sigma_e, sigma0,
                                                C = NULL, Z = NULL, W = NULL,
                                                var_Z = NULL, var_C = NULL, var_CT = NULL, mean_C = NULL,
                                                cov_WC = NULL, cov_WZ = NULL, mean_W = NULL, var_W = NULL, N = NULL, R = NULL,
                                                simple = TRUE){
  if(is.null(C) == TRUE){
    T <- length(var_CT)
    u <- sigma_e^2 + T * sigma0^2
    M1_11 <- T/u * (R^2 + R)/(R + 1)^3 * N
    M1_12 <- 1/(R + 1) * mean_C/u * N
    M1_22 <- 1/sigma_e^2 * (N - 1) * sum(var_CT) - sigma0^2/sigma_e^2/u * (var_C) * (N - 1)
    M1 <- matrix(c(M1_11, M1_12, M1_12, M1_22), nrow = 2)
    M <- M1
    if(simple == FALSE){
      w2 <- var_W * (N - 1)
      w1 <- mean_W
      lambda <- u / N^2 / T / w2
      v <- (N - 1)^2
      M2 <- lambda * v * c(N*T * cov_WZ/u, N*cov_WC/u) %*% t(c(N*T * cov_WZ/u, N*cov_WC/u))
      M <- M1 - M2
    }
  } else {
    T <- nrow(C)
    N <- ncol(C)
    ai <- calculate_ai(sigma_e = sigma_e, sigma0 = sigma0, sigma1 = 0, C)
    bi <- calculate_bi(sigma_e = sigma_e, sigma0 = sigma0, sigma1 = 0, C)
    u <- sigma_e^2 + T * sigma0^2

    M1_11 <- T/u * var(Z[1, ]) * (N - 1)
    M1_12 <- (N - 1) * cov(ai, Z[1, ])
    M1_22 <- 1/sigma_e^2 * sum((C - apply(C, 1, mean))^2) - sigma0^2/sigma_e^2/u * sum((ai - mean(ai))^2) * u^2
    M1 <- matrix(c(M1_11, M1_12, M1_12, M1_22), nrow = 2)
    M <- M1
    if(simple == FALSE){
      w2 <- var(W[1, ]) * (N - 1)
      w1 <- sum(W[1, ])
      lambda <- u / N^2 / T / w2
      v <- (N - 1)^2
      M2 <- lambda * (N - 1)^2 * c(N*T * cov(W[1, ]/u, Z[1, ]), N*cov(ai, W[1, ])) %*% t(c(N*T * cov(W[1, ]/u, Z[1, ]), N * cov(ai, W[1, ])))
      M <- M1 - M2
    }
  }
  return(M)
}


#' calculate_variance_random_slope
#'
#' This function is used to calculate the variance estimator of (theta_Z, theta_C) in random slope model.
#' There are several senario we want to calculate variance:
#' 1. Known C, Z, W: do not need to specific later parameters. Can choose simple to be TRUE or FALSE
#' 2. Known C, Z: do not need to specify later parameters. Can choose simple to be TRUE
#' 3. Without known C, Z, W: need to specify later variance and covariance parameters and also sample size N. Can choose simple to be TRUE without known information about W
#'
#'
#' @param sigma_e random error of residual
#' @param sigma_0 random intercept variability
#' @param sigma_1 random slope variability
#' @param C fixed effect associated with baseline
#' @param Z random intercept
#' @param W random slope
#' @param var_Z variance of treatment Z
#' @param var_Cstar variance of (sum of Concentration level/(1+gi)) of each Subject
#' @param var_CstarT vector of variance of (Concentration level/(1 + gi)) at each time point
#' @param mean_Cstar mean of (sum of Concentration level/(1 + gi)) of each subject
#' @param cov_WCstar covariance of (sum of Concentration level/(1 + gi)) and Baseline
#' @param cov_WZ covariance of treatment and Baseline
#' @param mean_W mean of baseline
#' @param var_W variance of baseline
#' @param N Number of subjects
#' @param R Ratio of treatment subjects v.s placebo subjects
#' @param simple TRUE: ignore correlation between treatment and baseline and between sum of Concentration level and Baseline
#' @return variance estimator of (theta_Z, theta_C)
#' @examples
#' @export
calculate_variance_random_slope <- function(sigma_e, sigma0, sigma1,
                                            C = NULL, Z = NULL, W = NULL, mean_gi2 = NULL,
                                            var_Z = NULL, var_Cstar = NULL, var_CstarT = NULL, mean_Cstar = NULL,
                                            cov_WA = NULL, cov_WCstar = NULL, mean_W = NULL, var_W = NULL, N = NULL, R = NULL,
                                            simple = TRUE, approx = FALSE){
  if(is.null(C) == TRUE){
    T <- length(var_CstarT)
    u <- sigma_e^2 + T * sigma0^2
    M1_11 <- T/u * (R^2 + R)/(R + 1)^3 * N
    M1_12 <- 1/(R + 1) * mean_Cstar/u * N
    if(sigma1 == 0){
      M1_22 <- 1/sigma_e^2 * (N - 1) * sum(var_CstarT) - sigma0^2/sigma_e^2/u * (var_Cstar) * (N - 1)
    } else{
      M1_22 <- 1/sigma1^2 * mean_gi2 * N +  1/sigma_e^2 * (N - 1) * sum(var_CstarT) - sigma0^2/sigma_e^2/u * (var_Cstar) * (N - 1)
    }
    M1 <- matrix(c(M1_11, M1_12, M1_12, M1_22), nrow = 2)
    M <- M1
    if(simple == FALSE){
      w2 <- var_W * (N - 1)
      w1 <- mean_W
      lambda <- u / N^2 / T / w2
      v <- (N - 1)^2
      M2 <- lambda * v * c(N*T * cov_WZ/u, N*cov_WCstar/u) %*% t(c(N*T * cov_WZ/u, N*cov_WCstar/u))
      M <- M1 - M2
    }
  } else {
      T <- nrow(C)
      N <- ncol(C)
      ai <- calculate_ai(sigma_e = sigma_e, sigma0 = sigma0, sigma1 = sigma1, C)
      bi <- calculate_bi(sigma_e = sigma_e, sigma0 = sigma0, sigma1 = sigma1, C)
      gi <- calculate_gi(sigma_e = sigma_e, sigma0 = sigma0, sigma1 = sigma1, C)
      u <- sigma_e^2 + T * sigma0^2
      w2 <- var(W[1, ]) * (N - 1)
      w1 <- sum(W[1, ])
      lambda <- u / N^2 / T / w2
      v <- (N - 1)^2

      M1_11 <- T/u * var(Z[1, ]) * (N - 1)
      M1_12 <- (N - 1) * cov(ai/(1 + gi), Z[1, ])
      # approximation
      if(approx == TRUE){
        M1_22 <- sum(1/sigma1^2 * gi^2/(1 + gi)^2)  + sum(1/sigma1^2 * gi/(1 + gi)^2) - (sigma_e^2/N * t(rowSums(t(t(bi)/(1 + gi)))) %*% rowSums(t(t(bi)/(1 + gi))) + sigma0^2/N * (sum(ai/(1 + gi)))^2)
        # ------------
        # this is old approximation using square root
        # M1_22 <- sum(1/sigma1^2 * gi/(1 + gi)) - (sigma_e^2/N * t(rowSums(t(t(bi)/sqrt(1 + gi)))) %*% rowSums(t(t(bi)/sqrt(1 + gi))) + sigma0^2/N * (sum(ai/sqrt(1 + gi)))^2)
        # ------------
      } else {
        M1_22 <- sum(1/sigma1^2 * gi/(1 + gi)) - (sigma_e^2/N * t(rowSums(t(t(bi)/(1 + gi)))) %*% rowSums(t(t(bi)/(1 + gi))) + sigma0^2/N * (sum(ai/(1 + gi)))^2) #true element
      }
      #(sum(apply(t(C)/sqrt(1 + gi), 2, sum)^2)/N/sigma_e^2 - sigma0^2/sigma_e^2/u/N * sum(colSums(t(C)/sqrt(1 + gi)))^2) - (sigma_e^2/N * t(rowSums(t(t(bi)/sqrt(1 + gi)))) %*% rowSums(t(t(bi)/sqrt(1 + gi))) + sigma0^2/N * (sum(ai/sqrt(1 + gi)))^2) check the approximation is valid

      M1 <- matrix(c(M1_11, M1_12, M1_12, M1_22), nrow = 2)
      M2 <- lambda * v * c(N*T * cov(W[1, ]/u, Z[1, ]), N*cov(ai/(1 + gi), W[1, ])) %*% t(c(N*T * cov(W[1, ]/u, Z[1, ]), N * cov(ai/(1 + gi), W[1, ])))
      if(simple == FALSE){
        M <- M1 - M2
      }else{
        M <- M1
      }
  }
  return(M)
}

