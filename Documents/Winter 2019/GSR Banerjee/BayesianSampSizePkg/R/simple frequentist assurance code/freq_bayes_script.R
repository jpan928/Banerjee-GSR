rm(list = ls())

# Script for Assurance Functions coded in Spring 2018
# good news is we're returning similar values as the original bayesian and frequentist function
# outputs


# Runs the Frequentist, Bayesian Assurance (in closed form), and saves them
# as two separate vectors before applying pwr_curve to plot the 2 sets of outcomes
# By the parameters we've selected, the results overlap perfectly


library(MASS)
library(ggplot2)
set.seed(1)

source("R/simple frequentist assurance code/pwr_curves.R")
source("R/simple frequentist assurance code/frequentist.R")
source("R/simple frequentist assurance code/bayes_assurance.R")

##########
# Case 1 #
##########
# Weak Priors: Overlapping Frequentist and Bayesian curves
sig_sq <- runif(1, 0, 1)
theta_0 <- K <- 0.15
theta_1 <- mu_beta_d <- runif(1,0,1)
n_d <- 1e+8
n_a <- 1e-8
a <- 1
Vbeta_d <- 1 / n_d # n_d should approach infinity, so Vbeta_d should be 0
Vbeta_a_inv <- 0 # n_a should approach 0
mu_beta_a <- 0
alpha <- 0.05

n <- seq(1, 200, 1)
y1 <- pwr_freq(n, sig_sq, alpha, theta_0, theta_1)
y2 <- assurance_nd_na(theta_0, theta_1, sig_sq, n, n_a, n_d, alpha)

png("R/images/curve1.png")
pwr_curve(n, y1, y2, TRUE, TRUE)
dev.off()



##########
# Case 2 #
##########
# Retain the strong design stage prior and adjust the analysis stage prior

# Case 2.1
n_a <- 10
n <- seq(1, 200, 1)
y1 <- pwr_freq(n, sig_sq, alpha, theta_0, theta_1)
y2 <- assurance_nd_na(theta_0, theta_1, sig_sq, n, n_a, n_d, alpha)

png("R/images/curve2.png")
pwr_curve(n, y1, y2, TRUE, TRUE)
dev.off()

# Case 2.2
n_a <- 100
n <- seq(1, 200, 1)
y1 <- pwr_freq(n, sig_sq, alpha, theta_0, theta_1)
y2 <- assurance_nd_na(theta_0, theta_1, sig_sq, n, n_a, n_d, alpha)

png("R/images/curve3.png")
pwr_curve(n, y1, y2, TRUE, TRUE)
dev.off()

# Case 2.3
n_a <- 500
n <- seq(1, 200, 1)
y1 <- pwr_freq(n, sig_sq, alpha, theta_0, theta_1)
y2 <- assurance_nd_na(theta_0, theta_1, sig_sq, n, n_a, n_d, alpha)

png("R/images/curve4.png")
pwr_curve(n, y1, y2, TRUE, TRUE)
dev.off()




##########
# Case 3 #
##########
# Remove the strong design prior and adjust the analysis stage prior
n_d <- 1
n_a <- 10
y1 <- pwr_freq(n, sig_sq, alpha, theta_0, theta_1)
y2 <- assurance_nd_na(theta_0, theta_1, sig_sq, n, n_a, n_d, alpha)

png("R/images/curve5.png")
pwr_curve(n, y1, y2, TRUE, TRUE)
dev.off()



##########
# Case 4 #
##########
# Weak Design Prior

n_d <- 1e-8

n <- seq(1, 200, 1)
y1 <- pwr_freq(n, sig_sq, alpha, theta_0, theta_1)
y2 <- assurance_nd_na(theta_0, theta_1, sig_sq, n, n_a, n_d, alpha)

png("R/images/curve6.png")
pwr_curve(n, y1, y2, TRUE, TRUE)
dev.off()










