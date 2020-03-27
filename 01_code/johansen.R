### global Parameters
df <- read_csv(here::here('00_data/lutkepohl.csv'))
df <- df%>%
    dplyr::select(linvestment,
                  lincome, 
                  lconsumption)


lags <- 1
trend <- 'const'

#### local Parameters (Ergebnisse gleich zu stata, daher wahrscheinlich auch in R lags + 1)

jlags <- lags + 1

#### test

#### just right hand variables 

jo_test <- urca::ca.jo(df, type = "eigen", K = jlags, ecdet = trend)

stat[1, 2] <- jo_test@lambda[1]
jo_nvar <- jo_test@P
AA <- jo_test@teststat
urca::summary(jo_test)



