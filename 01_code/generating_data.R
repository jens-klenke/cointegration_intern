##all cointehrated

## Simulation to test 

library("urca")

set.seed(123)

z <- rep(0, 10000)
for (i in 2:10000) z[i] <- z[i-1] + rnorm(1)

p <- q <- r <- rep(0, 10000)

p <- 0.3*z + rnorm(10000)
q <- 0.6*z + rnorm(10000)
r <- 0.2*z + rnorm(10000)

df <- data.frame(p,q,r)

#write.csv(df, file = '00_data/simulation.csv')

### just 2 integrated
set.seed(123456)

z <- rep(0, 10000)
for (i in 2:10000) z[i] <- z[i-1] + rnorm(1)

p <- q <- r <- rep(0, 10000)

p <- 0.3*z + rnorm(10000)
q <- 0.6*z + rnorm(10000)
r <-   rnorm(10000, 10, 100)

df_2 <- data.frame(p,q,r)


write.csv(df, file = '00_data/simulation_2.csv')
