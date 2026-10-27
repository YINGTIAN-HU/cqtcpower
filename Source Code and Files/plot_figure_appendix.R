# -------------------------------
# Appendix Figure Code
# -------------------------------

library(dplyr)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(cowplot)

# Figure B1
power <- read.table("./Results/power.txt")
power.correct <- power %>% filter((theta_Z_true ==0 & CV==0 & model == "RI") |
                                    (theta_Z_true ==0.5 & CV==0 & model == "RI + Z") |
                                    (theta_Z_true ==0 & CV==0.75 & model == "RS0") |
                                    (theta_Z_true ==0.5 & CV==0.75 & model == "RS0 + Z")
)
power.correct <- power.correct %>% filter(rho == 0)

power.correct$C_supra <- factor(power.correct$C_supra, levels = c(0.008, 0.0113), labels = c(TeX("$\\C_{supra} = 0.008$"),
                                                                                             TeX("$\\C_{supra} = 0.0113$")))
power.correct$eta0 <- factor(power.correct$eta0, levels = c(5, 8))
power.correct$epsilon <- factor(power.correct$epsilon, levels = c(5, 8)) 
power.correct$var <- interaction(power.correct$eta0, power.correct$epsilon)
levels(power.correct$var) <- c(TeX("$\\eta_{0} = 5, \\sigma_{\\epsilon} = 5$"),
                               TeX("$\\eta_{0} = 8,  \\sigma_{\\epsilon} = 5$"),
                               TeX("$\\eta_{0} = 5,  \\sigma_{\\epsilon} = 8$"),
                               TeX("$\\eta_{0} = 8,  \\sigma_{\\epsilon} = 8$"))

power.correct <- power.correct %>%
  pivot_longer(
    cols = c("analytical", "empiricalKR"),
    names_to = "type",
    values_to = "power") %>% 
  mutate(Method = factor(type, 
                         levels = c("analytical", "empiricalKR"),
                         labels = c("analytical", "empirical")))  

power.correct$Model <- power.correct$model
p <- ggplot(power.correct, aes(x = delta, y = power, group = interaction(Method, Model))) + 
  geom_line(aes(color = Method, linetype = Model)) +
  geom_point(aes(color = Method, shape = Model)) + 
  geom_abline(aes(intercept = 0.05, slope = 0)) +   
  facet_grid(C_supra ~ var, labeller=label_parsed) + xlab(expression(Delta*Delta*"QTc")) + 
  theme(
    strip.text = element_text(size=15, face = "bold"),
    axis.title.x = element_text(size=25, face = "bold"),
    axis.title.y = element_text(size=25, face = "bold"),
    axis.text.y = element_text(size = 15, face= "bold"),
    axis.text.x = element_text(size = 15, face= "bold"),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.spacing = unit(0.6, "lines"),
    plot.margin = margin(6, 6, 6, 6),
    legend.title = element_text(size=15, face= "bold"),
    legend.text = element_text(size=15)) + ylab("Power")
legend <- get_legend(
  # create some space to the left of the legend
  p + theme(legend.box.margin = margin(0, 0, 0, 12))
)


pdf("./Figure/power_comparison.pdf", width = 18, height = 8)
plot_grid(p + theme(legend.position = "none"), legend, rel_widths = c(4,.4))
dev.off()


# Figure B2
power.two <- power %>% filter((theta_Z_true ==0 & CV==0.75 & model == "RS0 + Z") |
                                (theta_Z_true ==0.5 & CV==0.75 & model == "RS + Z")
)
power.two <- power.two %>% filter(rho %in% c(0.2, 0.5, 0.8))

power.two$rho <- factor(power.two$rho, levels = c(0.2, 0.5, 0.8), labels = c(TeX("$\\rho = 0.2$"),
                                                                             TeX("$\\rho = 0.5$"),
                                                                             TeX("$\\rho = 0.8$")))
power.two$C_supra <- factor(power.two$C_supra, levels = c(0.008, 0.0113), labels = c(TeX("$\\C_{supra} = 0.008$"),
                                                                                     TeX("$\\C_{supra} = 0.0113$")))
power.two$eta0 <- factor(power.two$eta0, levels = c(5, 8))
power.two$epsilon <- factor(power.two$epsilon, levels = c(5, 8)) 
power.two$var <- interaction(power.two$eta0, power.two$epsilon)
levels(power.two$var) <- c(TeX("$\\eta_{0} = 5, \\sigma_{\\epsilon} = 5$"),
                           TeX("$\\eta_{0} = 8,  \\sigma_{\\epsilon} = 5$"),
                           TeX("$\\eta_{0} = 5,  \\sigma_{\\epsilon} = 8$"),
                           TeX("$\\eta_{0} = 8,  \\sigma_{\\epsilon} = 8$"))


power.two <- power.two %>%
  pivot_longer(
    cols = c("analytical", "empiricalKR"),
    names_to = "type",
    values_to = "Power") %>% 
  mutate(Method = factor(type, 
                         levels = c("analytical", "empiricalKR"),
                         labels = c("analytical", "empirical")))  

power.two$Model <- as.factor(power.two$model)
p <- ggplot(power.two, aes(x = delta, y = Power, group = interaction(Method, Model))) + 
  geom_line(aes(color = Method, linetype = Model)) +
  geom_point(aes(color = Method, shape = Model)) + 
  geom_abline(aes(intercept = 0.05, slope = 0)) +   
  facet_grid(C_supra + rho ~ var, labeller=label_parsed) + xlab(expression(Delta*Delta*"QTc")) + 
  theme(
    strip.text = element_text(size=18, face = "bold"),
    axis.title.x = element_text(size=25, face = "bold"),
    axis.title.y = element_text(size=25, face = "bold"),
    axis.text.y = element_text(size = 15, face= "bold"),
    axis.text.x = element_text(size = 15, face= "bold"),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.spacing = unit(0.6, "lines"),
    plot.margin = margin(6, 6, 6, 6),
    legend.title = element_text(size=15, face= "bold"),
    legend.text = element_text(size=15)) + ylab("Power")

legend <- get_legend(
  # create some space to the left of the legend
  p + theme(legend.box.margin = margin(0, 0, 0, 12))
)

pdf("./Figure/power_comparison2.pdf", width = 18, height = 24)
plot_grid(p + theme(legend.position = "none"), legend, rel_widths = c(4,.4))
dev.off()
