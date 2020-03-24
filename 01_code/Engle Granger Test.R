# global parameters 
lags <- 5

# Engle Granger Test
E_G_lm <- lm(p ~ q, data = df)
error <- E_G_lm$residuals
adf <- ur.df(error, type = "none", lags = lags)

stat[1, 1] <- adf@teststat[1]
