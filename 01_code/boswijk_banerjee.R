### Local Parameters
lags <- 4

# df <- data.frame(
#    x = c(11:20),
#    y = c(21:30))

# Xlag <- cbind(q, r)
Xlag <- cbind(depVar, indepVar)
Y_dif <- diff(depVar)
W <- diff(indepVar)

### Lagged values
Xlag_diff <- as.data.frame(diff(as.matrix(Xlag)))

x <- matrix(NA, nrow = nrow(Xlag_diff), ncol = lags)
X <- c()

if (lags >= 1) {
    for (j in 1:ncol(Xlag_diff)) {
        for (i in 1:lags) {
            x[, i] <- Hmisc::Lag(Xlag_diff[, j], shift = i)
        }
        X <- cbind(X, x)
    }
    W <- cbind(W, X)
}

W