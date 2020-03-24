### global Parameters
lags <- 2
trendtype <- 'trend'
#### local Parameters (Ergebnisse gleich zu stata, daher wahrscheinlich auch in R lags + 1)

jlags <- lags + 1

#### test

#### just right hand variables 

jo_test <- urca::ca.jo(df, type = "eigen", K = jlags, ecdet = trendtype)

stat[1, 2] <- jo_test@lambda[1]
jo_nvar <- jo_test@P

urca::summary(jo_test)

