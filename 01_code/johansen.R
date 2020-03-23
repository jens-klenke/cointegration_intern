### globals Parameters
lags <- 2
trendtype <- 'trend'
#### local Parameters (Ergenisse gleich zu stata, daher wahrscheilich auch in R lags + 1 )

jlags <- lags + 1

#### test

#### just right hand variables 

jo_test <- ca.jo(df, type="eigen", K = jlags, ecdet = trendtype)

stat[1,2] <- jo_test@lambda[1]
jo_nvar <- jo_test@P



