### Local Parameters
lags <- c(1, 2, 5)

Xlag <- cbind(q, r)
Y_dif <- diff(q)
W <- diff(r)

### Lagged values