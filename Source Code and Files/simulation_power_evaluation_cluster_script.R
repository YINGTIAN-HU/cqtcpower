# ----------------------------------------------
# cluster simulation script
# ----------------------------------------------
library(nlme)
library(lme4)
library(lmerTest)
library(MASS)

setwd("/projects/hu_lab/Yingtian/cQTC/Code")
source("./utils.R")
source("./generate_concentration.R")
source("./generate_delta_QTC.R")
source("./generate_parallel_SAD_data.R")
source("./calculate_variance.R")
source("./calculate_power.R")

args=(commandArgs(TRUE))
print(args)

if (length(args)==0) {
  U1=5
  CV=0
  theta_Z_true=0.5
  N_A=6
  N_P=2
  eta0=8
  epsilon=5
  rho = 0
  C_supra = 0.008
} else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

if(U1 == 10){
  N_sim <- 10000
} else {
  N_sim <- 100
}

res <- calculate_power(U1 = U1,
                       CV = CV, # CV of random slope
                       theta_Z_true = theta_Z_true,
                       C_supra = C_supra,
                       N_A = N_A,
                       N_P = N_P,
                       t = c(1, 2, 4, 8, 12, 16, 24),
                       dose_level = c(1, 3, 10, 30, 100, 200, 400, 800),
                       eta0 = eta0,
                       epsilon = epsilon,
                       rho = rho,
                       theta_T_true = c(1.93, 2.29, 3.53, -3.15, -0.28, 5.17, 3.98),
                       N_sim = N_sim)

folder_name <- paste("/projects/hu_lab/Yingtian/cQTC/Results/bias_U1", U1, "_CV", CV, "_theta_Z_true", theta_Z_true, "_N_A", N_A, "_N_P", N_P, "_eta0", eta0, "_epsilon", epsilon, "_rho", rho, "_C_supra", C_supra, sep = "")
command <- paste("mkdir -p ", folder_name, sep =  "")
system(command)
setwd(folder_name)
write.table(x = res$bias, file = "bias.txt")
write.table(x = res$CP, file = "CP.txt")
write.table(x = res$power, file = "power.txt")

