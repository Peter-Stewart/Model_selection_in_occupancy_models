# Set working directory to wherever the outputs of the simulations are stored ####
setwd("C:/temp/occupancy_paper")

# Load packages #####
library(MASS)
library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)

# Load results ####
# Scenario 1
resultsdf_m1_occ <- get(load("resultsdf_m1_occ.Rdata"))
resultsdf_m2_occ <- get(load("resultsdf_m2_occ.Rdata"))

predictions_all_occ <- get(load("predictions_all_occ.Rdata"))
retrodictions_all_occ <- get(load("retrodictions_all_occ.Rdata"))

Akaike_all_occ <- get(load("Akaike_all_occ.Rdata"))

# Scenario 2
resultsdf_m1_det <- get(load("resultsdf_m1_det.Rdata"))
resultsdf_m2_det <- get(load("resultsdf_m2_det.Rdata"))

predictions_all_det <- get(load("predictions_all_det.Rdata"))
retrodictions_all_det <- get(load("retrodictions_all_det.Rdata"))

Akaike_all_det <- get(load("Akaike_all_det.Rdata"))

# Scenario 3
resultsdf_m1_occ_det <- get(load("resultsdf_m1_occ_det.Rdata"))
resultsdf_m2_occ_det <- get(load("resultsdf_m2_occ_det.Rdata"))
resultsdf_m3_occ_det <- get(load("resultsdf_m3_occ_det.Rdata"))
resultsdf_m4_occ_det <- get(load("resultsdf_m4_occ_det.Rdata"))

predictions_all_occ_det <- get(load("predictions_all_occ_det.Rdata"))
retrodictions_all_occ_det <- get(load("retrodictions_all_occ_det.Rdata"))

Akaike_all_occ_det <- get(load("Akaike_all_occ_det.Rdata"))


# MAIN TEXT FIGURES ####
# Figure 4 ####
# Plots of true effect vs estimated effect, with colour showing whether 95% CI contains true value

# Scenario 1, model 1
sub1 <- resultsdf_m1_occ %>% select(B1, X.est, X.ci.low, X.ci.up, X.in95ci)
sub1_in <- sub1 %>% filter(X.in95ci==1)
sub1_out <- sub1 %>% filter(X.in95ci==0)
p1 <- ggplot(NULL, aes(x=B1, y=X.est)) +
  geom_point(data = sub1_out, shape=1, size=1.5) + 
  geom_point(data = sub1_in, shape=19, colour ="#426fb3", size=1.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(-1,1) + ylim(-10,10) +
  xlab(expression("True effect of X on "* psi)) + ylab(expression("Estimated effect of X on " * psi)) +
  ggtitle(expression(psi * " ~ " * X)) +
  theme_classic()

# Scenario 1, model 2
sub2 <- resultsdf_m2_occ %>% select(B1, X.est, X.ci.low, X.ci.up, X.in95ci)
sub2_in <- sub2 %>% filter(X.in95ci==1)
sub2_out <- sub2 %>% filter(X.in95ci==0)
p2 <- ggplot(NULL, aes(x=B1, y=X.est)) +
  geom_point(data = sub2_out, shape=1, size=1.5) + 
  geom_point(data = sub2_in, shape=19, colour ="#426fb3", size=1.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(-1,1) + ylim(-10,10) +
  xlab(expression("True effect of X on "* psi)) + ylab(expression("Estimated effect of X on " * psi)) +
  ggtitle(expression(psi * " ~ " * X + D)) +
  theme_classic()

# Scenario 2, model 1
sub3 <- resultsdf_m1_det %>% select(B1, X.est, X.ci.low, X.ci.up, X.in95ci)
sub3_in <- sub3 %>% filter(X.in95ci==1)
sub3_out <- sub3 %>% filter(X.in95ci==0)
p3 <- ggplot(NULL, aes(x=B1, y=X.est)) +
  geom_point(data = sub3_out, shape=1, size=1.5) + 
  geom_point(data = sub3_in, shape=19, colour ="#426fb3", size=1.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(-1,1) + ylim(-10,10) +
  xlab(expression("True effect of X on "* psi)) + ylab(expression("Estimated effect of X on " * psi)) +
  ggtitle(expression(psi * " ~ " * X * ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 2, model 2
sub4 <- resultsdf_m2_det %>% select(B1, X.est, X.ci.low, X.ci.up, X.in95ci)
sub4_in <- sub4 %>% filter(X.in95ci==1)
sub4_out <- sub4 %>% filter(X.in95ci==0)
p4 <- ggplot(NULL, aes(x=B1, y=X.est)) +
  geom_point(data = sub4_out, shape=1, size=1.5) + 
  geom_point(data = sub4_in, shape=19, colour ="#426fb3", size=1.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(-1,1) + ylim(-10,10) +
  xlab(expression("True effect of X on "* psi)) + ylab(expression("Estimated effect of X on " * psi)) +
  ggtitle(expression(psi * " ~ " * X * ", " * p * " ~ " * U + R)) +
  theme_classic()

# Scenario 3, model 1
sub5 <- resultsdf_m1_occ_det %>% select(B1, X.est, X.ci.low, X.ci.up, X.in95ci)
sub5_in <- sub5 %>% filter(X.in95ci==1)
sub5_out <- sub5 %>% filter(X.in95ci==0)
p5 <- ggplot(NULL, aes(x=B1, y=X.est)) +
  geom_point(data = sub5_out, shape=1, size=1.5) + 
  geom_point(data = sub5_in, shape=19, colour ="#426fb3", size=1.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(-1,1) + ylim(-10,10) +
  xlab(expression("True effect of X on "* psi)) + ylab(expression("Estimated effect of X on " * psi)) +
  ggtitle(expression(psi * " ~ " * X * ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 3, model 2
sub6 <- resultsdf_m2_occ_det %>% select(B1, X.est, X.ci.low, X.ci.up, X.in95ci)
sub6_in <- sub6 %>% filter(X.in95ci==1)
sub6_out <- sub6 %>% filter(X.in95ci==0)
p6 <- ggplot(NULL, aes(x=B1, y=X.est)) +
  geom_point(data = sub6_out, shape=1, size=1.5) + 
  geom_point(data = sub6_in, shape=19, colour ="#426fb3", size=1.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(-1,1) + ylim(-10,10) +
  xlab(expression("True effect of X on "* psi)) + ylab(expression("Estimated effect of X on " * psi)) +
  ggtitle(expression(psi * " ~ " * X + D* ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 3, model 3
sub7 <- resultsdf_m3_occ_det %>% select(B1, X.est, X.ci.low, X.ci.up, X.in95ci)
sub7_in <- sub7 %>% filter(X.in95ci==1)
sub7_out <- sub7 %>% filter(X.in95ci==0)
p7 <- ggplot(NULL, aes(x=B1, y=X.est)) +
  geom_point(data = sub7_out, shape=1, size=1.5) + 
  geom_point(data = sub7_in, shape=19, colour ="#426fb3", size=1.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(-1,1) + ylim(-10,10) +
  xlab(expression("True effect of X on "* psi)) + ylab(expression("Estimated effect of X on " * psi)) +
  ggtitle(expression(psi * " ~ " * X * ", " * p * " ~ " * U + R)) +
  theme_classic()

# Scenario 3, model 4
sub8 <- resultsdf_m4_occ_det %>% select(B1, X.est, X.ci.low, X.ci.up, X.in95ci)
sub8_in <- sub8 %>% filter(X.in95ci==1)
sub8_out <- sub8 %>% filter(X.in95ci==0)
p8 <- ggplot(NULL, aes(x=B1, y=X.est)) +
  geom_point(data = sub8_out, shape=1, size=1.5) + 
  geom_point(data = sub8_in, shape=19, colour ="#426fb3", size=1.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(-1,1) + ylim(-10,10) +
  xlab(expression("True effect of X on "* psi)) + ylab(expression("Estimated effect of X on " * psi)) +
  ggtitle(expression(psi * " ~ " * X + D* ", " * p * " ~ " * U + R)) +
  theme_classic()

# Arrange plots in grid
fig4 <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8,
                  nrow=2,
                  labels = c("A)", "B)", "C)", "D)", "E)", "F)", "G)", "H)"))
fig4

# Figure 5 ####
# Kernel density contours for predictive accuracy (mean absolute error vs. proportion of sites in 95% CI)

# Scenario 1, model 1 
p9 <- ggplot(predictions_all_occ, aes(x=m1_prop_in95ci, y=m1_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) +
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X)) +
  theme_classic()

# Scenario 1, model 2
p10 <- ggplot(predictions_all_occ, aes(x=m2_prop_in95ci, y=m2_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X + D)) +
  theme_classic()

# Scenario 2, model 1
# Need to set bandwidth manually as results are all clustered in one region
bandwidth.nrd(predictions_all_det$m1_mean_abs_error)
bandwidth.nrd(predictions_all_det$m1_prop_in95ci)
p11 <- ggplot(predictions_all_det, aes(x=m1_prop_in95ci, y=m1_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE, h=c(0.001, 0.005163218)) +
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 2, model 2
# Need to set bandwidth manually as results are all clustered in one region
bandwidth.nrd(predictions_all_det$m2_mean_abs_error)
bandwidth.nrd(predictions_all_det$m2_prop_in95ci)
p12 <- ggplot(predictions_all_det, aes(x=m2_prop_in95ci, y=m2_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE, h=c(0.001, 0.005090337)) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U + R)) +
  theme_classic()

# Scenario 3, model 1
p13 <- ggplot(predictions_all_occ_det, aes(x=m1_prop_in95ci, y=m1_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 3, model 2
p14 <- ggplot(predictions_all_occ_det, aes(x=m2_prop_in95ci, y=m2_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X + D* ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 3, model 3
p15 <- ggplot(predictions_all_occ_det, aes(x=m3_prop_in95ci, y=m3_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U + R)) +
  theme_classic()

# Scenario 3, model 4
p16 <- ggplot(predictions_all_occ_det, aes(x=m4_prop_in95ci, y=m4_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X + D* ", " * p * " ~ " * U + R)) +
  theme_classic()

# Arrange plots in grid
fig5 <- plot_grid(p9, p10, p11, p12, p13, p14, p15, p16,
                  nrow=2,
                  labels = c("A)", "B)", "C)", "D)", "E)", "F)", "G)", "H)"))
fig5


# Figure 6 ####
# AIC and BIC weights for each model
# First need to order simulations by weight
ord1 <- Akaike_all_occ[order(Akaike_all_occ$m1_weight),]
ord1b <- Akaike_all_occ[order(Akaike_all_occ$m1_BIC_weight),]
ord2 <- Akaike_all_occ[order(Akaike_all_occ$m2_weight),]
ord2b <- Akaike_all_occ[order(Akaike_all_occ$m2_BIC_weight),]

ord3 <- Akaike_all_det[order(Akaike_all_det$m1_weight),]
ord3b <- Akaike_all_det[order(Akaike_all_det$m1_BIC_weight),]
ord4 <- Akaike_all_det[order(Akaike_all_det$m2_weight),]
ord4b <- Akaike_all_det[order(Akaike_all_det$m2_BIC_weight),]

ord5 <- Akaike_all_occ_det[order(Akaike_all_occ_det$m1_weight),]
ord5b <- Akaike_all_occ_det[order(Akaike_all_occ_det$m1_BIC_weight),]
ord6 <- Akaike_all_occ_det[order(Akaike_all_occ_det$m2_weight),]
ord6b <- Akaike_all_occ_det[order(Akaike_all_occ_det$m2_BIC_weight),]
ord7 <- Akaike_all_occ_det[order(Akaike_all_occ_det$m3_weight),]
ord7b <- Akaike_all_occ_det[order(Akaike_all_occ_det$m3_BIC_weight),]
ord8 <- Akaike_all_occ_det[order(Akaike_all_occ_det$m4_weight),]
ord8b <- Akaike_all_occ_det[order(Akaike_all_occ_det$m4_BIC_weight),]

# Scenario 1, model 1
p17 <- ggplot(NULL, aes(x=1:1000)) +
  ylim(0,1) +
  geom_area(aes(y=ord1b$m1_BIC_weight), alpha=1, fill="#24878e") +
  geom_area(aes(y=ord1$m1_weight), alpha=1, fill="#fde725") +
  geom_line(aes(y=ord1$m1_weight), size=0.8, colour="black") +
  geom_line(aes(y=ord1b$m1_BIC_weight), size=0.8, colour="black") +
  ylab("Weight") + xlab("Simulation") +
  geom_abline(slope=0, intercept=0, lty=2, size=1) + 
  geom_abline(slope=0, intercept=1, lty=2, size=1) +
  geom_abline(slope=0, intercept=0.5, lty=3, size=1) +
  ggtitle(expression(psi * " ~ " * X)) +
  theme_classic()

# Scenario 1, model 2
p18 <- ggplot(NULL, aes(x=1:1000)) +
  ylim(0,1) +
  geom_area(aes(y=ord2$m2_weight), alpha=1, fill="#fde725") +
  geom_area(aes(y=ord2b$m2_BIC_weight), alpha=1, fill="#24878e") +
  geom_line(aes(y=ord2$m2_weight), size=0.8, colour="black") +
  geom_line(aes(y=ord2b$m2_BIC_weight), size=0.8, colour="black") +
  ylab("Weight") + xlab("Simulation") +
  geom_abline(slope=0, intercept=0, lty=2, size=1) + 
  geom_abline(slope=0, intercept=1, lty=2, size=1) +
  geom_abline(slope=0, intercept=0.5, lty=3, size=1) +
  ggtitle(expression(psi * " ~ " * X + D)) +
  theme_classic()


# Scenario 2, model 1
p19 <- ggplot(NULL, aes(x=1:1000)) +
  ylim(0,1) +
  geom_area(aes(y=ord3b$m1_BIC_weight), alpha=1, fill="#24878e") +
  geom_area(aes(y=ord3$m1_weight), alpha=1, fill="#fde725") +
  geom_line(aes(y=ord3$m1_weight), size=0.8, colour="black") +
  geom_line(aes(y=ord3b$m1_BIC_weight), size=0.8, colour="black") +
  ylab("Weight") + xlab("Simulation") +
  geom_abline(slope=0, intercept=0, lty=2, size=1) + 
  geom_abline(slope=0, intercept=1, lty=2, size=1) +
  geom_abline(slope=0, intercept=0.5, lty=3, size=1) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 2, model 2
p20 <- ggplot(NULL, aes(x=1:1000)) +
  ylim(0,1) +
  geom_area(aes(y=ord4$m2_weight), alpha=1, fill="#fde725") +
  geom_area(aes(y=ord4b$m2_BIC_weight), alpha=1, fill="#24878e") +
  geom_line(aes(y=ord4$m2_weight), size=0.8, colour="black") +
  geom_line(aes(y=ord4b$m2_BIC_weight), size=0.8, colour="black") +
  ylab("Weight") + xlab("Simulation") +
  geom_abline(slope=0, intercept=0, lty=2, size=1) + 
  geom_abline(slope=0, intercept=1, lty=2, size=1) +
  geom_abline(slope=0, intercept=0.5, lty=3, size=1) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U + R)) +
  theme_classic()

# Scenario 3, model 1
p21 <- ggplot(NULL, aes(x=1:1000)) +
  ylim(0,1) +
  geom_area(aes(y=ord5$m1_weight), alpha=1, fill="#fde725") +
  geom_area(aes(y=ord5b$m1_BIC_weight), alpha=1, fill="#24878e") +
  geom_line(aes(y=ord5$m1_weight), size=0.8, colour="black") +
  geom_line(aes(y=ord5b$m1_BIC_weight), size=0.8, colour="black") +
  ylab("Weight") + xlab("Simulation") +
  geom_abline(slope=0, intercept=0, lty=2, size=1) + 
  geom_abline(slope=0, intercept=1, lty=2, size=1) +
  geom_abline(slope=0, intercept=0.5, lty=3, size=1) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 3, model 2
p22 <- ggplot(NULL, aes(x=1:1000)) +
  ylim(0,1) +
  geom_area(aes(y=ord6b$m2_BIC_weight), alpha=1, fill="#24878e") +
  geom_area(aes(y=ord6$m2_weight), alpha=1, fill="#fde725") +
  geom_line(aes(y=ord6$m2_weight), size=0.8, colour="black") +
  geom_line(aes(y=ord6b$m2_BIC_weight), size=0.8, colour="black") +
  ylab("Weight") + xlab("Simulation") +
  geom_abline(slope=0, intercept=0, lty=2, size=1) + 
  geom_abline(slope=0, intercept=1, lty=2, size=1) +
  geom_abline(slope=0, intercept=0.5, lty=3, size=1) +
  ggtitle(expression(psi * " ~ " * X + D* ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 3, model 3
p23 <- ggplot(NULL, aes(x=1:1000)) +
  ylim(0,1) +
  geom_area(aes(y=ord7b$m3_BIC_weight), alpha=1, fill="#24878e") +
  geom_area(aes(y=ord7$m3_weight), alpha=1, fill="#fde725") +
  geom_line(aes(y=ord7$m3_weight), size=0.8, colour="black") +
  geom_line(aes(y=ord7b$m3_BIC_weight), size=0.8, colour="black") +
  ylab("Weight") + xlab("Simulation") +
  geom_abline(slope=0, intercept=0, lty=2, size=1) + 
  geom_abline(slope=0, intercept=1, lty=2, size=1) +
  geom_abline(slope=0, intercept=0.5, lty=3, size=1) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U + R)) +
  theme_classic()

# Scenario 3, model 4
p24 <- ggplot(NULL, aes(x=1:1000)) +
  ylim(0,1) +
  geom_area(aes(y=ord8$m4_weight), alpha=1, fill="#fde725") +
  geom_area(aes(y=ord8b$m4_BIC_weight), alpha=1, fill="#24878e") +
  geom_line(aes(y=ord8$m4_weight), size=0.8, colour="black") +
  geom_line(aes(y=ord8b$m4_BIC_weight), size=0.8, colour="black") +
  ylab("Weight") + xlab("Simulation") +
  geom_abline(slope=0, intercept=0, lty=2, size=1) + 
  geom_abline(slope=0, intercept=1, lty=2, size=1) +
  geom_abline(slope=0, intercept=0.5, lty=3, size=1) +
  ggtitle(expression(psi * " ~ " * X + D* ", " * p * " ~ " * U + R)) +
  theme_classic()

# Arrange plots in grid
fig6 <- plot_grid(p17, p18, p19, p20, p21, p22, p23, p24,
                  nrow=2,
                  labels = c("A)", "B)", "C)", "D)", "E)", "F)", "G)", "H)"))
fig6


# SUPPLEMENTARY FIGURES ####
# Figure S1 ####
# Kernel density contours for retrodictive accuracy (mean absolute error vs. proportion of sites in 95% CI)
# Scenario 1, model 1
ps1a <- ggplot(retrodictions_all_occ, aes(x=m1_prop_in95ci, y=m1_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) +
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X)) +
  theme_classic()

# Scenario 1, model 2
ps1b <- ggplot(retrodictions_all_occ, aes(x=m2_prop_in95ci, y=m2_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X + D)) +
  theme_classic()

# Scenario 2, model 1
# Need to set bandwidth manually as all simulations are clustered in one area of the graph
bandwidth.nrd(retrodictions_all_det$m1_mean_abs_error)
bandwidth.nrd(retrodictions_all_det$m1_prop_in95ci)
ps1c <- ggplot(retrodictions_all_det, aes(x=m1_prop_in95ci, y=m1_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE, h=c(0.001, 0.005165063)) +
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 2, model 2
# Need to set bandwidth manually as all simulations are clustered in one area of the graph
bandwidth.nrd(retrodictions_all_det$m2_mean_abs_error)
bandwidth.nrd(retrodictions_all_det$m2_prop_in95ci)
ps1d <- ggplot(retrodictions_all_det, aes(x=m2_prop_in95ci, y=m2_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE, h=c(0.001, 0.005091841)) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U + R)) +
  theme_classic()

# Scenario 3, model 1
ps1e <- ggplot(retrodictions_all_occ_det, aes(x=m1_prop_in95ci, y=m1_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 3, model 2
ps1f <- ggplot(retrodictions_all_occ_det, aes(x=m2_prop_in95ci, y=m2_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X + D* ", " * p * " ~ " * U)) +
  theme_classic()

# Scenario 3, model 3
ps1g <- ggplot(retrodictions_all_occ_det, aes(x=m3_prop_in95ci, y=m3_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X* ", " * p * " ~ " * U + R)) +
  theme_classic()

# Scenario 3, model 4
ps1h <- ggplot(retrodictions_all_occ_det, aes(x=m4_prop_in95ci, y=m4_mean_abs_error)) +
  geom_density2d_filled(show.legend = FALSE) + 
  ylab("Mean absolute error") + xlab("Proportion in 95% C.I.") +
  xlim(0, 1) + ylim(0, 0.2) +
  ggtitle(expression(psi * " ~ " * X + D* ", " * p * " ~ " * U + R)) +
  theme_classic()

# Arrange plots in grid
figs1 <- plot_grid(ps1a, ps1b, ps1c, ps1d, ps1e, ps1f, ps1g, ps1h, 
                   nrow=2,
                   labels = c("A)", "B)", "C)", "D)", "E)", "F)", "G)", "H)"))
figs1


# Figure S2 ####
# Pairs plot showing relationship between effect sizes and AIC/BIC weights for scenario 1, model 1
b1 <- resultsdf_m1_occ$B1
b2 <- resultsdf_m1_occ$B2
b3 <- resultsdf_m1_occ$B3
b4 <- resultsdf_m1_occ$B4
b5 <- resultsdf_m1_occ$B5
check <- cbind(b1, b2, b3, b4, b5, Akaike_all_occ)
check2 <- check %>% select(b1, b2, b3, b4, b5, m1_weight, m1_BIC_weight)
pairs(check2, upper.panel = NULL, labels = c(expression(beta*X*psi),
                                             expression(beta*A*X),
                                             expression(beta*A*D),
                                             expression(beta*C*D),
                                             expression(beta*C*psi)
                                             ,"Akaike weight", "BIC weight"))

# Figure S3 ####
# Histograms of parameter values when estimate of X effect is in 95% CI for scenario 1, model 2
check <- resultsdf_m2_occ %>% select(B1, B2, B3, B4, B5, X.in95ci) %>% filter(X.in95ci==1)
par(mfrow=c(2,3))
hist(check$B1, breaks = 50, main = expression(beta*X*psi), xlab = "Value") 
hist(check$B2, breaks = 50, main = expression(beta*A*X), xlab = "Value") 
hist(check$B3, breaks = 50, main = expression(beta*A*D), xlab = "Value") 
hist(check$B4, breaks = 50, main = expression(beta*C*D), xlab = "Value")
hist(check$B5, breaks = 50, main = expression(beta*C*psi), xlab = "Value")
par(mfrow=c(1,1))

# Figure S4 ####
# Histograms of parameter values when estimate of X effect is in 95% CI for scenario 3, model 2
check <- resultsdf_m2_occ_det %>% select(B1, B2, B3, B4, B5, X.in95ci) %>% filter(X.in95ci==1)
par(mfrow=c(2,3))
hist(check$B1, breaks = 50, main = expression(beta*X*psi), xlab = "Value") 
hist(check$B2, breaks = 50, main = expression(beta*A*X), xlab = "Value") 
hist(check$B3, breaks = 50, main = expression(beta*A*D), xlab = "Value") 
hist(check$B4, breaks = 50, main = expression(beta*C*D), xlab = "Value")
hist(check$B5, breaks = 50, main = expression(beta*C*psi), xlab = "Value")
par(mfrow=c(1,1))

# Figure S5 ####
# Histograms of parameter values when estimate of X effect is in 95% CI for scenario 3, model 4
check <- resultsdf_m4_occ_det %>% select(B1, B2, B3, B4, B5, X.in95ci) %>% filter(X.in95ci==1)
par(mfrow=c(2,3))
hist(check$B1, breaks = 50, main = expression(beta*X*psi), xlab = "Value") 
hist(check$B2, breaks = 50, main = expression(beta*A*X), xlab = "Value") 
hist(check$B3, breaks = 50, main = expression(beta*A*D), xlab = "Value") 
hist(check$B4, breaks = 50, main = expression(beta*C*D), xlab = "Value")
hist(check$B5, breaks = 50, main = expression(beta*C*psi), xlab = "Value")
par(mfrow=c(1,1))

# Figure S6 ####
# Pairs plot showing relationship between effect sizes and AIC/BIC weights for scenario 3, model 3
b1 <- resultsdf_m3_occ_det$B1
b2 <- resultsdf_m3_occ_det$B2
b3 <- resultsdf_m3_occ_det$B3
b4 <- resultsdf_m3_occ_det$B4
b5 <- resultsdf_m3_occ_det$B5
check <- cbind(b1, b2, b3, b4, b5, Akaike_all_occ_det)
check2 <- check %>% select(b1, b2, b3, b4, b5, m3_weight, m3_BIC_weight)
pairs(check2, upper.panel = NULL, labels = c(expression(beta*X*psi),
                                             expression(beta*A*X),
                                             expression(beta*A*D),
                                             expression(beta*C*D),
                                             expression(beta*C*psi)
                                             ,"Akaike weight", "BIC weight"))
