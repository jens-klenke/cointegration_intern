# stat 

stat <- matrix(NA, nrow = 1, ncol =  4)

# global parameters 
lags <- 5

# Engle Granger Test
E_G_lm <- lm(p ~ q, data = df)
error <- E_G_lm$residuals
adf <- urca::ur.df(error, type = "none", lags = lags)

stat[1, 1] <- adf@teststat[1]