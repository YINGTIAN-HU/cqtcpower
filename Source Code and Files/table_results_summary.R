# ---------------------
# Summarize Simluation Results; Table 4, 5, 6 results
# ---------------------
CV <- 0
N_A <- 6
N_P <- 2
eta0 <- 8
epsilon <- 5
rho <- 0
C_supra <- 0.008

U1 <- 5
theta_Z_true <- 0

summary.results <- function(CV = 0 , N_A = 6, N_P = 2, eta0 = 9, epsilon = 5, rho = 0, C_supra = 0.008, U1 = 5, theta_Z_true = 0){
  folder.name <- "./Results/"
  subfolder.name <- paste("bias_U1", U1, "_CV", CV, "_theta_Z_true", theta_Z_true, "_N_A", N_A, "_N_P", N_P, "_eta0", eta0, "_epsilon", epsilon, "_rho", rho, "_C_supra", C_supra, "/",sep = "")
  folder.full.name <- paste0(folder.name, subfolder.name)
  bias <- round(read.table(paste0(folder.full.name, "bias.txt")) * 100, 3)
  CP <- round(read.table(paste0(folder.full.name, "CP.txt")), 3)
  power <- round(read.table(paste0(folder.full.name, "power.txt")), 3)
  dat <- cbind(power[, c(1,2)], CP[,1], bias)
  colnames(dat) <- c("analytical", "empirical", "CP", "bias")
  return(dat)
}

# ------
# Table 4
# ------
summary.results(CV = 0 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0, C_supra = 0.008, U1 = 10, theta_Z_true = 0)
summary.results(CV = 0 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0, C_supra = 0.008, U1 = 5, theta_Z_true = 0)
summary.results(CV = 0 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0, C_supra = 0.008, U1 = 10, theta_Z_true = 0.5)
summary.results(CV = 0 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0, C_supra = 0.008, U1 = 5, theta_Z_true = 0.5)

# ------
# Table 5
# ------
summary.results(CV = 0.75 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0, C_supra = 0.008, U1 = 10, theta_Z_true = 0.5)
summary.results(CV = 0.75 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0, C_supra = 0.008, U1 = 5, theta_Z_true = 0.5)
summary.results(CV = 0.25 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0, C_supra = 0.008, U1 = 10, theta_Z_true = 0.5)
summary.results(CV = 0.25 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0, C_supra = 0.008, U1 = 5, theta_Z_true = 0.5)

# ------
# Table 6
# ------
summary.results(CV = 0.75 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0.2, C_supra = 0.008, U1 = 10, theta_Z_true = 0.5)
summary.results(CV = 0.75 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0.5, C_supra = 0.008, U1 = 10, theta_Z_true = 0.5)
summary.results(CV = 0.75 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0.8, C_supra = 0.008, U1 = 10, theta_Z_true = 0.5)
summary.results(CV = 0.75 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0.2, C_supra = 0.008, U1 = 5, theta_Z_true = 0.5)
summary.results(CV = 0.75 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0.5, C_supra = 0.008, U1 = 5, theta_Z_true = 0.5)
summary.results(CV = 0.75 , N_A = 6, N_P = 2, eta0 = 8, epsilon = 5, rho = 0.8, C_supra = 0.008, U1 = 5, theta_Z_true = 0.5)
