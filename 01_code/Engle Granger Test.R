library(urca)
simulation_2 <- read.csv("simulation_2.csv")

# global parameters 
lags <- 2

# Engle Granger Test
E_G_lm <- lm(p ~ q, data = simulation_2)
error <- E_G_lm$residuals
adf <- ur.df(error, type = "none", lags = lags)

library(aTSA)
coint.test