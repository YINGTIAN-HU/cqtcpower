library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

source("utils.R")
source("calculate_analytical_sample_size.R")
source("calculate_variance.R")

delta_seq <- c(10, 9, 8, 7, 6, 5)
#baseline
S0 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N_G = 8, #number of subject in each group
                                           R = 3, #N_A: N_P
                                           Cmax = 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01),
                                           delta = delta_seq[i],
                                           sigma_e  = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S0[i, ] <- unlist(res)
}

#S1
S1 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N_G = 8, #number of subject in each group
                                           R = 3, #N_A: N_P
                                           Cmax = 1.5 * 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01),
                                           delta = delta_seq[i],
                                           sigma_e  = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S1[i, ] <- unlist(res)
}
#S2
S2 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N_G = 8, #number of subject in each group
                                           R = 3, #N_A: N_P
                                           Cmax = 0.75 * 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01),
                                           delta = delta_seq[i],
                                           sigma_e  = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S2[i, ] <- unlist(res)
}

#S3
S3 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N_G = 8, #number of subject in each group
                                           R = 3, #N_A: N_P
                                           Cmax = 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01),
                                           delta = delta_seq[i],
                                           sigma_e  = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S3[i, ] <- unlist(res)
}

#S4
S4 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N_G = 8, #number of subject in each group
                                           R = 3, #N_A: N_P
                                           Cmax = 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01),
                                           delta = delta_seq[i],
                                           sigma_e = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S4[i, ] <- unlist(res)
}

#S5
S5 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N_G = 8, #number of subject in each group
                                           R = 3, #N_A: N_P
                                           Cmax = 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.45, 0.8, 0.6, 0.3, 0.15, 0.05, 0.01, 0.005),
                                           delta = delta_seq[i],
                                           sigma_e = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S5[i, ] <- unlist(res)
}

#S6
S6 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N_G = 8, #number of subject in each group
                                           R = 3, #N_A: N_P
                                           Cmax = 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.45, 0.9, 0.95, 0.8, 0.4, 0.05),
                                           delta = delta_seq[i],
                                           sigma_e  = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S6[i, ] <- unlist(res)
}

#S7
S7 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N_G = 8, #number of subject in each group
                                           R = 3, #N_A: N_P
                                           Cmax = 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.9, 0.95, 0.85, 0.4),
                                           delta = delta_seq[i],
                                           sigma_e  = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S7[i, ] <- unlist(res)
}

######################################
# #S8
# S8 <- array(0, dim = c(4, 4))
# for(i in 1:length(delta_seq)){
#   res <- calculate_analytical_sample_size2(N = 8, #number of subject in each group
#                            R = 3, #N_A: N_P
#                            Cmax = 0.008, #maximum averaged concetration level
#                            C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
#                            profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1),
#                            delta = delta_seq[i],
#                            sigma_e  = 6,
#                            sigma0 = 3,
#                            C_supra = 0.008,
#                            theta0 = 0)
#   S8[i, ] <- unlist(res)
# }
#
# #S9
# S9 <- rep(0, 4)
# for(i in 1:length(delta_seq)){
#   res <- calculate_analytical_sample_size2(N = 8, #number of subject in each group
#                            R = 3, #N_A: N_P
#                            Cmax = 0.008, #maximum averaged concetration level
#                            C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
#                            profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4),
#                            delta = delta_seq[i],
#                            sigma_e  = 6,
#                            sigma0 = 3,
#                            C_supra = 0.008,
#                            theta0 = 0)
#   S9[i] <- res[2]
# }
#
# #S10
# S10 <- rep(0, 4)
# for(i in 1:length(delta_seq)){
#   res <- calculate_analytical_sample_size2(N = 8, #number of subject in each group
#                            R = 3, #N_A: N_P
#                            Cmax = 0.008, #maximum averaged concetration level
#                            C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
#                            profile = c(0.45, 0.9, 0.95, 0.85, 0.75),
#                            delta = delta_seq[i],
#                            sigma_e  = 6,
#                            sigma0 = 3,
#                            C_supra = 0.008,
#                            theta0 = 0)
#   S10[i] <- res[2]
# }
#
# #S11
# S11 <- rep(0, 4)
# for(i in 1:length(delta_seq)){
#   res <- calculate_analytical_sample_size2(N = 8, #number of subject in each group
#                            R = 3, #N_A: N_P
#                            Cmax = 0.008, #maximum averaged concetration level
#                            C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
#                            profile = c(0.45, 0.9, 0.95, 0.85),
#                            delta = delta_seq[i],
#                            sigma_e  = 6,
#                            sigma0 = 3,
#                            C_supra = 0.008,
#                            theta0 = 0)
#   S11[i] <- res[2]
# }
#
# #S12
# S12 <- rep(0, 4)
# for(i in 1:length(delta_seq)){
#   res <- calculate_analytical_sample_size2(N = 8, #number of subject in each group
#                            R = 3, #N_A: N_P
#                            Cmax = 0.008, #maximum averaged concetration level
#                            C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
#                            profile = c(0.45, 0.9, 0.95),
#                            delta = delta_seq[i],
#                            sigma_e  = 6,
#                            sigma0 = 3,
#                            C_supra = 0.008,
#                            theta0 = 0)
#   S12[i] <- res[2]
# }
#
#
# #S13
# S13 <- rep(0, 4)
# for(i in 1:length(delta_seq)){
#   res <- calculate_analytical_sample_size2(N = 8, #number of subject in each group
#                            R = 3, #N_A: N_P
#                            Cmax = 0.008, #maximum averaged concetration level
#                            C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
#                            profile = c(0.45, 0.9),
#                            delta = delta_seq[i],
#                            sigma_e  = 6,
#                            sigma0 = 3,
#                            C_supra = 0.008,
#                            theta0 = 0)
#   S13[i] <- res[2]
# }
#
# #S14
# S14 <- rep(0, 4)
# for(i in 1:length(delta_seq)){
#   res <- calculate_analytical_sample_size2(N = 8, #number of subject in each group
#                            R = 3, #N_A: N_P
#                            Cmax = 0.008, #maximum averaged concetration level
#                            C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
#                            profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.1, 0.01),
#                            delta = delta_seq[i],
#                            sigma_e  = 6,
#                            sigma0 = 3,
#                            C_supra = 0.008,
#                            theta0 = 0)
#   S14[i] <- res[2]
# }
#########################################

S <- data.frame("S0" =  matrix(S0, ncol =1, byrow = TRUE))/8
S$model <- rep(c("RI", "RI + Z", "RS0", "RS0 + Z"), each = 6)
S$delta <- rep(10:5, 4)
S$S1 <- matrix(S1, ncol =1, byrow = TRUE)/8
S$S2 <- matrix(S2, ncol =1, byrow = TRUE)/8
S$S3 <- matrix(S5, ncol =1, byrow = TRUE)/8
S$S4 <- matrix(S6, ncol =1, byrow = TRUE)/8
S$S5 <- matrix(S7, ncol =1, byrow = TRUE)/8
S <- S[, c(2, 3, 1, 4:8)]
S <- S %>%
  pivot_longer(
    cols = starts_with("S"),
    names_to = "Scenario",
    values_to = "Size")

S$delta <- 10 - S$delta
S$delta <- factor(S$delta)

plot_model <- function(model_name){
  ggplot(S %>% dplyr::filter(model == model_name), aes(x = delta, y = Size, fill = Scenario)) + geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual("Scenario", values = c("S0" = "#F8766D",
                                             "S1" = "orange",
                                             "S2" = "#00B81F",
                                             "S3" = "cyan",
                                             "S4" = "#00A5FF",
                                             "S5" = "#E76BF3")) + ylim(c(0, 30)) + scale_x_discrete(name = expression(Delta*Delta*"QTC"),
                                                                                                    limits=factor(c(0, 1, 2, 3, 4, 5))) +
    theme(
      strip.text = element_text(size=25),
      axis.title.x = element_text(size=35, ),
      axis.title.y = element_text(size=35),
      axis.text.y = element_text(size = 15, face= "bold"),
      axis.text.x = element_text(size = 15, face= "bold"),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      plot.margin = margin(6, 6, 6, 6),
      legend.title = element_text(size=13, face= "bold"),
      legend.text = element_text(size=13)) #+ geom_text(x="1", y=250, label=model_name, size = 5)
  
}

p1 <- plot_model("RI") + ylab(NULL)
p2 <- plot_model("RI + Z") + ylab(NULL)
p3 <- plot_model("RS0") + ylab(NULL)
p4 <- plot_model("RS0 + Z") + ylab(NULL)
# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# # arrange the three plots in a single row
# prow <- plot_grid(
#   p1 + theme(legend.position="none", axis.title.x=element_blank()) + ggtitle("RI") + xlab(NULL),
#   p2 + theme(legend.position="none") + ggtitle("RI + Z")  + ggtitle("RI + Z") + xlab(NULL),
#   p3 + theme(legend.position="none") + ggtitle("RS0")  + ggtitle("RS0") + xlab(NULL),
#   p4 + theme(legend.position="none") + ggtitle("RS0 + Z")  + ggtitle("RS0 + Z") + xlab(NULL),
#   #labels = c("RI", "RI + Z", "RS0", "RS0 + Z"),
#   #hjust = -2,
#   #vjust = 2,
#   nrow = 2
# )
# prow

yleft <- textGrob("Sample Size Per Cohort",
                  rot = 90, gp = gpar(fontsize = 20))
bottom <- textGrob(expression(Delta*Delta*"QTc"), gp = gpar(fontsize = 20))

p1 <- p1 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RI")
p2 <- p2 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RI + Z")  + ggtitle("RI + Z")
p3 <- p3 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RS0")  + ggtitle("RS0") + xlab(NULL)
p4 <- p4 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RS0 + Z")  + ggtitle("RS0 + Z")

#p = list(p1,p2,p3,p4) %>% map(~.x + labs(x=NULL, y=NULL))
uni <- grid.arrange(grobs=list(p1, p2, p3, p4), ncol = 4, nrow = 1,
                    left = yleft, bottom = bottom)



# add the legend to the row we made earlier. Give it one-third of
# the width of one plot (via rel_widths).
pdf("./Figure/Figure 2.pdf", width = 18, height = 5)
plot_grid(uni, legend, rel_widths = c(4, .22))
dev.off()



S <- data.frame("S0" =  matrix(S0, ncol =1, byrow = TRUE))
S$model <- rep(c("RI", "RI + Z", "RS0", "RS0 + Z"), each = 6)
S$delta <- rep(10:5, 4)
S$S6 <- matrix(S3, ncol =1, byrow = TRUE)
S$S7 <- matrix(S4, ncol =1, byrow = TRUE)
S <- S[, c(2, 1, 3:5)]
S <- S %>%
  pivot_longer(
    cols = starts_with("S"),
    names_to = "Scenario",
    values_to = "Size")

S$delta <- 10 - S$delta
S$delta <- factor(S$delta)

plot_model <- function(model_name){
  ggplot(S %>% dplyr::filter(model == model_name), aes(x = delta, y = Size, fill = Scenario)) + geom_bar(stat="identity", position=position_dodge()) +
    #scale_fill_manual("Scenario", values = c("S0" = "#F8766D",
    #                                        "S6" = "orange",
    #                                        "S7" = "#00B81F")) +
    ylim(c(0, 200)) + scale_x_discrete(name = expression(Delta*Delta*"QTc"), limits=factor(c(0, 1, 2, 3, 4, 5))) +
    theme(
      strip.text = element_text(size=25),
      axis.title.x = element_text(size=35),
      axis.title.y = element_text(size=35),
      axis.text.y = element_text(size = 15, face= "bold"),
      axis.text.x = element_text(size = 15, face= "bold"),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      plot.margin = margin(6, 6, 6, 6),
      legend.title = element_text(size=13, face= "bold"),
      legend.text = element_text(size=13)) #+ geom_text(x="1", y=250, label=model_name, size = 5)
  
}

p1 <- plot_model("RI") + ylab(NULL)
p2 <- plot_model("RI + Z") + ylab(NULL)
p3 <- plot_model("RS0") + ylab(NULL)
p4 <- plot_model("RS0 + Z") + ylab(NULL)
# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

yleft <- textGrob("Total Sample Size",
                  rot = 90, gp = gpar(fontsize = 20))
bottom <- textGrob(expression(Delta*Delta*"QTc"), gp = gpar(fontsize = 20))

p1 <- p1 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RI")
p2 <- p2 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RI + Z")  + ggtitle("RI + Z")
p3 <- p3 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RS0")  + ggtitle("RS0") + xlab(NULL)
p4 <- p4 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RS0 + Z")  + ggtitle("RS0 + Z")

#p = list(p1,p2,p3,p4) %>% map(~.x + labs(x=NULL, y=NULL))
uni <- grid.arrange(grobs=list(p1, p2, p3, p4), ncol = 4, nrow = 1,
                    left = yleft, bottom = bottom)


S <- data.frame("S0" =  matrix(S0, ncol =1, byrow = TRUE))/8
S$model <- rep(c("RI", "RI + Z", "RS0", "RS0 + Z"), each = 6)
S$delta <- rep(10:5, 4)
S$S6 <- matrix(S3, ncol =1, byrow = TRUE)/6
S$S7 <- matrix(S4, ncol =1, byrow = TRUE)/4
S <- S[, c(2, 1, 3:5)]
S <- S %>%
  pivot_longer(
    cols = starts_with("S"),
    names_to = "Scenario",
    values_to = "Size")

S$delta <- 10 - S$delta
S$delta <- factor(S$delta)

plot_model <- function(model_name){
  ggplot(S %>% dplyr::filter(model == model_name), aes(x = delta, y = Size, fill = Scenario)) + geom_bar(stat="identity", position=position_dodge()) +
    #scale_fill_manual("Scenario", values = c("S0" = "#F8766D",
    #                                        "S6" = "orange",
    #                                        "S7" = "#00B81F")) +
    ylim(c(0, 40)) + scale_x_discrete(name = expression(Delta*Delta*"QTc"), limits=factor(c(0, 1, 2, 3, 4, 5))) +
    theme(
      strip.text = element_text(size=25),
      axis.title.x = element_text(size=35, ),
      axis.title.y = element_text(size=35),
      axis.text.y = element_text(size = 15, face= "bold"),
      axis.text.x = element_text(size = 15, face= "bold"),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      plot.margin = margin(6, 6, 6, 6),
      legend.title = element_text(size=13, face= "bold"),
      legend.text = element_text(size=13)) #+ geom_text(x="1", y=250, label=model_name, size = 5)
  
}


p1 <- plot_model("RI") + ylab(NULL)
p2 <- plot_model("RI + Z") + ylab(NULL)
p3 <- plot_model("RS0") + ylab(NULL)
p4 <- plot_model("RS0 + Z") + ylab(NULL)
# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

yleft <- textGrob("Sample Size Per Cohort",
                  rot = 90, gp = gpar(fontsize = 20))
bottom <- textGrob(expression(Delta*Delta*"QTc"), gp = gpar(fontsize = 20))

p1 <- p1 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RI")
p2 <- p2 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RI + Z")  + ggtitle("RI + Z")
p3 <- p3 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RS0")  + ggtitle("RS0") + xlab(NULL)
p4 <- p4 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RS0 + Z")  + ggtitle("RS0 + Z")

#p = list(p1,p2,p3,p4) %>% map(~.x + labs(x=NULL, y=NULL))
uni2 <- grid.arrange(grobs=list(p1, p2, p3, p4), ncol = 4, nrow = 1,
                     left = yleft, bottom = bottom)

# add the legend to the row we made earlier. Give it one-third of
# the width of one plot (via rel_widths).
pdf("./Figure/Figure 3.pdf", width = 18, height = 10)
plot_grid(plot_grid(uni2, uni, nrow = 2, labels = c("(A)", "(B)")), legend, rel_widths = c(4,.22))
dev.off()




S8 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N = 6, #number of subject in each group
                                           R = 2, #N_A: N_P
                                           Cmax = 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01),
                                           delta = delta_seq[i],
                                           sigma_e  = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S8[i, ] <- unlist(res)
}

#S9
S9 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N = 5, #number of subject in each group
                                           R = 4, #N_A: N_P
                                           Cmax = 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01),
                                           delta = delta_seq[i],
                                           sigma_e  = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S9[i, ] <- unlist(res)
}

#S10
S10 <- array(0, dim = c(length(delta_seq), 4))
for(i in 1:length(delta_seq)){
  res <- calculate_analytical_sample_size2(N = 6, #number of subject in each group
                                           R = 5, #N_A: N_P
                                           Cmax = 0.008, #maximum averaged concetration level
                                           C_alpha = c(1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1), #proportion of maximum Cmax in each group
                                           profile = c(0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01),
                                           delta = delta_seq[i],
                                           sigma_e  = 5,
                                           sigma0 = 8,
                                           C_supra = 0.008,
                                           theta0 = 0)
  S10[i, ] <- unlist(res)
}

S <- data.frame("S0" =  matrix(S0, ncol =1, byrow = TRUE))/8
S$model <- rep(c("RI", "RI + Z", "RS0", "RS0 + Z"), each = 6)
S$delta <- rep(10:5, 4)
S$S8 <- matrix(S8, ncol =1, byrow = TRUE)/8
S$S9 <- matrix(S9, ncol =1, byrow = TRUE)/8
S$S10 <- matrix(S10, ncol =1, byrow = TRUE)/8

S <- S[, c(2, 3, 1, 4:6)]
S <- S %>%
  pivot_longer(
    cols = starts_with("S"),
    names_to = "Scenario",
    values_to = "Size")
colnames(S)[4] <- "Size"

S$delta <- 10 - S$delta
S$delta <- factor(S$delta)
S$Scenario <- factor(S$Scenario, levels = c("S8", "S0", "S9", "S10"))

plot_model <- function(model_name){
  ggplot(S %>% dplyr::filter(model == model_name), aes(x = delta, y = Size, fill = Scenario)) + geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual("Scenario", values = c("S0" = "#F8766D",
                                             "S8" = "orange",
                                             "S9" = "#00BA38",
                                             "S10" = "#619CFF")) + ylim(c(0, 25)) + scale_x_discrete(name = expression(Delta*Delta*"QTC"),
                                                                                                     limits=factor(c(0, 1, 2, 3, 4, 5))) +
    theme(
      strip.text = element_text(size=25),
      axis.title.x = element_text(size=35, face= "bold"),
      axis.title.y = element_text(size=35, face= "bold"),
      axis.text.y = element_text(size = 15, face= "bold"),
      axis.text.x = element_text(size = 15, face= "bold"),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      plot.margin = margin(6, 6, 6, 6),
      legend.title = element_text(size=13, face= "bold"),
      legend.text = element_text(size=13)) #+ geom_text(x="1", y=250, label=model_name, size = 5)
  
}

p1 <- plot_model("RI") + ylab(NULL)
p2 <- plot_model("RI + Z") + ylab(NULL)
p3 <- plot_model("RS0") + ylab(NULL)
p4 <- plot_model("RS0 + Z") + ylab(NULL)
# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

yleft <- textGrob("Sample Size Per Cohort",
                  rot = 90, gp = gpar(fontsize = 20))
bottom <- textGrob(expression(Delta*Delta*"QTc"), gp = gpar(fontsize = 20))

p1 <- p1 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RI")
p2 <- p2 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RI + Z")  + ggtitle("RI + Z")
p3 <- p3 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RS0")  + ggtitle("RS0") + xlab(NULL)
p4 <- p4 + theme(legend.position="none", axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle("RS0 + Z")  + ggtitle("RS0 + Z")

#p = list(p1,p2,p3,p4) %>% map(~.x + labs(x=NULL, y=NULL))
uni <- grid.arrange(grobs=list(p1, p2, p3, p4), ncol = 4, nrow = 1,
                    left = yleft, bottom = bottom)

pdf("./Figure/Figure 4.pdf", width = 18, height = 5)
plot_grid(uni, legend, rel_widths = c(4, 0.22))
dev.off()






















