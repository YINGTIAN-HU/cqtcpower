# -----------------------
# function used to summarize power across different simluation scenarios
# -----------------------

delta_seq <- c(3, 5, 7, 10)
CV_seq <- c(0, 0.25, 0.75)
theta_Z_true_seq <- c(0, 0.5)
N_A_seq <- 6
N_P_seq <- 2
eta0_seq <- c(5, 8)
epsilon_seq <- c(5, 8)
rho_seq <- c(0, 0.2, 0.5, 0.8)
C_supra_seq <- c(0.008, 0.0113)

power <- NULL
CP <- NULL
bias <- NULL
for(delta in delta_seq){
  for(CV in CV_seq){
    for(theta_Z_true in theta_Z_true_seq){
      for(N_A in N_A_seq){
        for(N_P in N_P_seq){
          for(eta0 in eta0_seq){
            for(epsilon in epsilon_seq){
              for(rho in rho_seq){
                for(C_supra in C_supra_seq){
                  filename <- paste("./Results/0728/bias_U1",delta,"_CV",CV,"_theta_Z_true", theta_Z_true,
                                    "_N_A", N_A, "_N_P", N_P, "_eta0", eta0, "_epsilon", epsilon, "_rho", rho, "_C_supra", C_supra, "/power.txt", sep = "")
                  
                  if( file.exists(filename)){
                    power.temp <- read.table(filename)
                    power.temp$delta <- delta
                    power.temp$CV <- CV
                    power.temp$theta_Z_true <- theta_Z_true
                    power.temp$N_A <- N_A
                    power.temp$N_P <- N_P
                    power.temp$eta0 <- eta0
                    power.temp$epsilon <- epsilon
                    power.temp$rho <- rho
                    power.temp$C_supra <- C_supra
                    power.temp$model <- c("RI", "RI + Z", "RS0", "RS0 + Z", "RS", "RS + Z")
                    power <- rbind(power, power.temp)
                  }
                  
                  filename <- paste("./Results/0728/bias_U1",delta,"_CV",CV,"_theta_Z_true", theta_Z_true,
                                    "_N_A", N_A, "_N_P", N_P, "_eta0", eta0, "_epsilon", epsilon, "_rho", rho, "_C_supra", C_supra, "/CP.txt", sep = "")
                  
                  if( file.exists(filename)){
                    CP.temp <- read.table(filename)
                    CP.temp <- read.table(filename)
                    CP.temp$delta <- delta
                    CP.temp$CV <- CV
                    CP.temp$theta_Z_true <- theta_Z_true
                    CP.temp$N_A <- N_A
                    CP.temp$N_P <- N_P
                    CP.temp$eta0 <- eta0
                    CP.temp$epsilon <- epsilon
                    CP.temp$rho <- rho
                    CP.temp$C_supra <- C_supra
                    CP.temp$model <- c("RI", "RI + Z", "RS0", "RS0 + Z", "RS", "RS + Z")
                    CP <- rbind(CP, CP.temp)
                  }
                  
                  filename <- paste("./Results/0728/bias_U1",delta,"_CV",CV,"_theta_Z_true", theta_Z_true,
                                    "_N_A", N_A, "_N_P", N_P, "_eta0", eta0, "_epsilon", epsilon, "_rho", rho, "_C_supra", C_supra, "/bias.txt", sep = "")
                  if( file.exists(filename)){
                    bias.temp <- read.table(filename)
                    bias.temp <- read.table(filename)
                    bias.temp <- read.table(filename)
                    bias.temp$delta <- delta
                    bias.temp$CV <- CV
                    bias.temp$theta_Z_true <- theta_Z_true
                    bias.temp$N_A <- N_A
                    bias.temp$N_P <- N_P
                    bias.temp$eta0 <- eta0
                    bias.temp$epsilon <- epsilon
                    bias.temp$rho <- rho
                    bias.temp$C_supra <- C_supra
                    bias.temp$model <- c("RI", "RI + Z", "RS0", "RS0 + Z", "RS", "RS + Z")
                    bias <- rbind(bias, bias.temp)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
power <- data.frame(power)
write.table(power, file = "./Results/power.txt")