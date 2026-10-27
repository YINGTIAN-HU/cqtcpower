library(nlme)
library(lme4)
library(lmerTest)
library(MASS)
library(scales)
library(ggplot2)
library(dplyr)
library(cowplot)

source("./generate_parallel_SAD_data.R")
source("./generate_delta_QTC.R")
source("./generate_concentration.R")

set.seed(0)
res <- parallel_SAD_sim(theta0 = -5.76, #
                                  theta_C = 9.5/0.008,
                                  theta_Z = 0.5,
                                  theta_W = -0.58,
                                  theta_T = c(1.93, 2.29, 3.53, -3.15, -0.28, 5.17, 3.98),
                                  eta0 = 8,
                                  eta1 = 0.75 * 9.5/0.008,
                                  epsilon = 5,
                                  mean_baseline = 408,
                                  sd_baseline = 16,
                                  N_A = 6,
                                  N_P = 2,
                                  t = c(1, 2, 4, 8, 12, 16, 24),
                                  dose_level = c(1, 3, 10, 30, 100, 200, 400, 800),
                                  C_supra = 0.008,
                                  fit_model = FALSE,
                                  model = "slope",
                                  rho = 0, sd = 0.1, method = 1)
res$data[1:5, 1:5]
data <- res$data
N <- length(unique(data$Subject))
Time <- nrow(data)/N
data$Time <- rep(c(1, 2, 4, 8, 12, 16, 24), N)
colnames(data) <- toupper(colnames(data))
data <- data %>% mutate(DOSE = DOSE * TREATMENT,
                        DOSE = as.factor(DOSE),
                        TIME = as.numeric(TIME),
                        TREATMENT = as.factor(TREATMENT))

data$CONCENTRATION <- data$CONCENTRATION * 1000
pdf(file = "./Figure/QTC.pdf", width = 20, height = 6.5)
p1 <- ggplot(data = data,
       aes(x = TIME, y = CONCENTRATION, group = SUBJECT, color = DOSE)) + geom_point() + geom_line() + facet_wrap(. ~ DOSE)  +
           scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                         labels = trans_format("log10", math_format(10^.x))) +
       ylab("CONCENTRATION (ug/mL)") + xlab("Time (h)") + theme_bw() + theme(legend.position = "none",
                                                                             strip.text = element_text(size=15),
                                                                             axis.title.x = element_text(size=18),
                                                                             axis.title.y = element_text(size=18),
                                                                             axis.text.y = element_text(size = 12, face= "bold"),
                                                                             axis.text.x = element_text(size = 12, face= "bold"),
                                                                             axis.line = element_line(colour = "black"),
                                                                             panel.grid.major = element_blank(),
                                                                             panel.grid.minor = element_blank(),
                                                                             panel.background = element_blank(),
                                                                             panel.spacing = unit(0.6, "lines"))




p2 <- ggplot(data = data, aes(x = TIME, y = DELTA_QTC, group = SUBJECT, color = DOSE)) + geom_point() + geom_line() +
      facet_wrap(. ~ DOSE)  + ylab(expression(Delta*"QTC")) + xlab("Time (h)") + theme_bw() + theme(legend.position = "none",
                                                                                                    strip.text = element_text(size=15),
                                                                                                    axis.title.x = element_text(size=18),
                                                                                                    axis.title.y = element_text(size=18),
                                                                                                    axis.text.y = element_text(size = 12, face= "bold"),
                                                                                                    axis.text.x = element_text(size = 12, face= "bold"),
                                                                                                    axis.line = element_line(colour = "black"),
                                                                                                    panel.grid.major = element_blank(),
                                                                                                    panel.grid.minor = element_blank(),
                                                                                                    panel.background = element_blank(),
                                                                                                    panel.spacing = unit(0.6, "lines"))
p3 <- ggplot(data = data, aes(x = CONCENTRATION, y = DELTA_QTC, group = TREATMENT, color= TREATMENT)) + geom_point() +
      geom_smooth(method='lm', color = "black", linetype = 2)  +
      ylab(expression(Delta*"QTC")) + xlab("CONCENTRATION (ug/mL)") + theme_bw() + theme(legend.position = "none",
                                                                                         strip.text = element_text(size=15),
                                                                                         axis.title.x = element_text(size=18, ),
                                                                                         axis.title.y = element_text(size=18),
                                                                                         axis.text.y = element_text(size = 12, face= "bold"),
                                                                                         axis.text.x = element_text(size = 12, face= "bold"),
                                                                                         axis.line = element_line(colour = "black"),
                                                                                         panel.grid.major = element_blank(),
                                                                                         panel.grid.minor = element_blank(),
                                                                                         panel.background = element_blank(),
                                                                                         panel.spacing = unit(0.6, "lines"))

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 20, nrow = 1)
dev.off()


# T <- 7
# C <- matrix(res$data$Concentration, nrow = T)
# W <- matrix(res$data$Baseline, nrow = T)
# Z <- matrix(res$data$Treatment, nrow = T)
# 
# # true variance estimator
# sigma_e <- 5 # sd of random error
# sigma0 <- 8 # sd of random intercept
# sigma1 <- 0 # sd of random slope
# calculate_variance_random_intercept(sigma_e  = sigma_e, sigma0 = sigma0, C = C, Z = Z, W = W, simple = TRUE)
# 
# # fit LME model
# formula <- as.formula(paste("delta_QTC ~ -1 + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "), " + Treatment + Concentration",  sep = ""))
# fit <- lme(formula, data = res$data, random = ~ 1|Subject, control = lmeControl(opt = "optim"), method = "ML")  #random intecept
# variance <- vcov(fit)[c(9,10), c(9,10)]
# variance
# 
# # variance estimator based derived formula
# sigma_e <- as.numeric(VarCorr(fit)[2, 2]) # sd of random error
# sigma0 <- as.numeric(VarCorr(fit)[1, 2]) # sd of random intercept
# sigma1 <- 0 # sd of random slope
# var_est <- calculate_variance_random_intercept(sigma_e  = sigma_e, sigma0 = sigma0, C = C, Z = Z, W = W, simple = FALSE)
# solve(var_est) - variance
# (solve(var_est) - variance)/variance
# 
# 
# # fit LME model
# formula <- as.formula(paste("delta_QTC ~ -1 + Baseline + ", paste("t", 1:T, sep = "", collapse = " + "), " + Treatment + Concentration - 1 + (Concentration||Subject)",  sep = ""))
# fit2 <- lmer(formula,
#             data = res$data, REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
# variance <- vcov(fit2)[c(9,10), c(9,10)]
# variance
# 
# # variance estimator based derived formula
# sigma_e <- attr(VarCorr(fit2), "sc") # sd of random error
# sigma0 <- sqrt(as.numeric(VarCorr(fit2)[1])) # sd of random intercept
# sigma1 <- sqrt(as.numeric(VarCorr(fit2)[2])) # sd of random slope
# var_est <- calculate_variance_random_slope(sigma_e  = sigma_e, sigma0 = sigma0, sigma1 = sigma1, C = C, Z = Z, W = W, simple = FALSE, approx = FALSE)
# var_est2 <- calculate_variance_random_slope(sigma_e  = sigma_e, sigma0 = sigma0, sigma1 = sigma1, C = C, Z = Z, W = W, simple = FALSE, approx = TRUE)
# 
# 
# solve(var_est) - variance
# (solve(var_est) - variance)/variance
# 
# solve(var_est2) - variance
# (solve(var_est2) - variance)/variance 

