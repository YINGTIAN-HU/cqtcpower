#' Calculation of gi
#'
#' This function is used to calculate gi, which is used when calculating the inverse of covariance matrix
#'
#' @param sigma_e random error of residual
#' @param sigma_0 random intercept variability
#' @param sigma_1 random slope variability
#' @param C fixed effect associated with baseline
#' @return variance estimator of (theta_Z, theta_C)
#' @examples
#' var <- calculate_variance_random_intercept(sigma_e, sigma0, C, Z, W, simple = TRUE)
#' var
#' @export
calculate_gi <- function(sigma_e, sigma0, sigma1, C){
  T <- nrow(C)
  N <- ncol(C)
  u <- sigma_e^2 + T * sigma0^2
  gi <- rep(0, N)
  for(i in 1:N){
    gi[i] <- sigma1^2/sigma_e^2 * sum(C[,i]^2) - sigma0^2 * sigma1^2/(sigma_e^2 * u) * (sum(C[,i]))^2
  }
  return(gi)
}


#' Calculation of ai
#'
#' This function is used to calculate ai, defined as \sum_kC_{ik}/(\sigma_e^2 + T\sigma_0^2)
#'
#'
#' @param sigma_e random error of residual
#' @param sigma_0 random intercept variability
#' @param sigma_1 random slope variability
#' @param C fixed effect associated with baseline
#' @return variance estimator of (theta_Z, theta_C)
#' @examples
#' var <- calculate_variance_random_intercept(sigma_e, sigma0, C, Z, W, simple = TRUE)
#' var
#' @export
calculate_ai <- function(sigma_e, sigma0, sigma1, C){
  T <- nrow(C)
  N <- ncol(C)
  u <- sigma_e^2 + T * sigma0^2
  ai <- rep(0, N)
  for(i in 1:N){
    ai[i] <- 1/u * sum(C[, i])
  }
  return(ai)
}

#' Calculation of bi
#'
#' This function is used to calculate ai, defined as \sum_kC_{ik}/(\sigma_e^2 + T\sigma_0^2)
#'
#' @param sigma_e random error of residual
#' @param sigma_0 random intercept variability
#' @param sigma_1 random slope variability
#' @param C fixed effect associated with baseline
#' @return variance estimator of (theta_Z, theta_C)
#' @examples
#' var <- calculate_variance_random_intercept(sigma_e, sigma0, C, Z, W, simple = TRUE)
#' var
#' @export
calculate_bi <- function(sigma_e, sigma0, sigma1, C){
  T <- nrow(C)
  N <- ncol(C)
  J <- rep(1, T)
  u <- sigma_e^2 + T * sigma0^2
  bi <- array(0, dim = c(T, N))
  for(i in 1:N){
    bi[, i] <- 1/sigma_e^2 * C[,i] - sigma0^2/(sigma_e^2 * u) * (sum(C[,i])) * J
  }
  return(bi)
}


calculate_CVC <- function(sigma_e, sigma0, sigma1, C, A, B){
  gi <- calculate_gi(sigma_e, sigma0, sigma1, C)
  res <- sum(1/sigma1^2 * gi/(1 + gi))
  return(res)
}
calculate_ZVZ <- function(sigma_e, sigma0, sigma1, C, A, B, Z){
  T <- nrow(C)
  N <- ncol(C)
  ai <- calculate_ai(sigma_e, sigma0, sigma1, C)
  gi <- calculate_gi(sigma_e, sigma0, sigma1, C)
  u <- sigma_e^2 + T * sigma0^2
  res <- sum((T/u - sigma1^2/(1 + gi)*ai^2) * Z[1, ]^2)
  return(res)
}
calculate_ZVE <- function(sigma_e, sigma0, sigma1, C, A, B, Z, W){
  T <- nrow(C)
  N <- ncol(C)
  J <- rep(1, T)
  ai <- calculate_ai(sigma_e, sigma0, sigma1, C)
  gi <- calculate_gi(sigma_e, sigma0, sigma1, C)
  bi <- calculate_bi(sigma_e, sigma0, sigma1, C)
  u <- sigma_e^2 + T * sigma0^2
  res <- c(sum((T/u - sigma1^2/(1 + gi)*ai^2) * Z[1, ] *W[1, ]),
           sum(1/u * Z[1, ]) * J - colSums(sigma1^2/(1 + gi)*ai*Z[1, ] * t(bi))
  )
  return(res)
}
calculate_ZVC <- function(sigma_e, sigma0, sigma1, C, A, B, Z, W){
  T <- nrow(C)
  N <- ncol(C)
  ai <- calculate_ai(sigma_e, sigma0, sigma1, C)
  gi <- calculate_gi(sigma_e, sigma0, sigma1, C)
  bi <- calculate_bi(sigma_e, sigma0, sigma1, C)
  u <- sigma_e^2 + T * sigma0^2
  res <- sum(1/(1+gi) * Z[1,] * ai)
  return(res)
}
calculate_EVC <- function(sigma_e, sigma0, sigma1, C, A, B, Z, W){
  T <- nrow(C)
  N <- ncol(C)
  ai <- calculate_ai(sigma_e, sigma0, sigma1, C)
  gi <- calculate_gi(sigma_e, sigma0, sigma1, C)
  bi <- calculate_bi(sigma_e, sigma0, sigma1, C)
  u <- sigma_e^2 + T * sigma0^2
  res <- c( sum(1/(1+gi) * ai *W[1,]), colSums(1/(1+gi) * t(bi)))
  return(res)
}

calculate_EVE <- function(sigma_e, sigma0, sigma1, C, A, B, Z, W){
  M11 <- calculate_M11(sigma_e, sigma0, sigma1, C, A, B, Z, W)
  M21 <- calculate_M21(sigma_e, sigma0, sigma1, C, A, B, Z, W)
  M22 <- calculate_M22(sigma_e, sigma0, sigma1, C, A, B, Z, W)
  M <- rbind(M11, M21)
  M2 <- rbind(t(M21), M22)
  M <- cbind(M, M2)
  return(M)
}

calculate_EVE_inv <- function(sigma_e, sigma0, sigma1, C, A, B, Z, W){
  T <- nrow(C)
  N <- ncol(C)
  J <- rep(1, T)
  JJ <- J %*% t(J)
  u <- sigma_e^2 + T * sigma0^2
  w2 <- var(W[1, ]) * (N - 1)
  w1 <- sum(W[1, ])

  alpha <- u / N^2 / T / w2
  M11 <-  N^2
  M12 <- -N * w1 * t(J)
  M21 <- -N * w1 * J
  M22 <- N*T*sigma_e^2/u * w2 * diag(T)  + N*T*sigma0^2/u * w2 * JJ + w1^2 * JJ
  M1 <- c(M11, M12)
  M2 <- cbind(M21, M22)
  M <- alpha * rbind(M1, M2)
  return(M)
}

calculate_M22 <- function(sigma_e, sigma0, sigma1, C, A, B, Z, W){
  T <- nrow(C)
  N <- ncol(C)
  ai <- calculate_ai(sigma_e, sigma0, sigma1, C)
  gi <- calculate_gi(sigma_e, sigma0, sigma1, C)
  bi <- calculate_bi(sigma_e, sigma0, sigma1, C)
  u <- sigma_e^2 + T * sigma0^2
  A_inv <- calculate_A_inv(sigma_e, sigma0, C)
  res <- array(0, dim = c(T, T))
  for(i in 1:N){
    #print(sigma1^2/(1 + gi[i]) * bi[, i] %*% t(bi[, i]))
    res <- res + A_inv - sigma1^2/(1 + gi[i]) * bi[, i] %*% t(bi[, i])
  }
  return(res)
}

calculate_M21 <- function(sigma_e, sigma0, sigma1, C, A, B, Z, W){
  T <- nrow(C)
  N <- ncol(C)
  ai <- calculate_ai(sigma_e, sigma0, sigma1, C)
  gi <- calculate_gi(sigma_e, sigma0, sigma1, C)
  bi <- calculate_bi(sigma_e, sigma0, sigma1, C)
  u <- sigma_e^2 + T * sigma0^2
  res <- array(0, dim = c(T, 1))
  for(i in 1:N){
    temp <- 1/u -  sigma1^2/(1 + gi[i]) * 1 * ai[i] * bi[, i]
    res <- res + W[1,i] * temp
  }
  return(res)
}

calculate_M11 <- function(sigma_e, sigma0, sigma1, C, A, B, Z, W){
  T <- nrow(C)
  N <- ncol(C)
  ai <- calculate_ai(sigma_e, sigma0, sigma1, C)
  gi <- calculate_gi(sigma_e, sigma0, sigma1, C)
  bi <- calculate_bi(sigma_e, sigma0, sigma1, C)
  u <- sigma_e^2 + T * sigma0^2
  res <- 0
  for(i in 1:N){
    res <- res + (T/u - sigma1^2/(1 + gi[i])*ai[i]^2) * t(W[1, i]) %*% W[1, i]
  }
  return(res)
}

calculate_A_inv <- function(sigma_e, sigma0, C){
  T <- nrow(C)
  N <- ncol(C)
  J <- rep(1, T)
  A_inv <- 1/sigma_e^2 * diag(T) - sigma0^2/(sigma_e^2 * (sigma_e^2 + T * sigma0^2)) * J %*% t(J)
  return(A_inv)
}
