### Local Parameters
lags <- c(1, 2, 5)

# df <- data.frame(
#    x = c(11:20),
#    y = c(21:30))

# Xlag <- cbind(q, r)
Xlag <- cbind(depVar, indepVar)
Y_dif <- diff(depVar) # Problem: 0 Variablen, muss wohl numeric werden (lapply?)
W <- diff(indepVar) # Problem: 0 Variablen, muss wohl numeric werden (lapply?)

### Lagged values
Xlag_diff <- as.data.frame(diff(as.matrix(Xlag)))

x <- matrix(NA, nrow = nrow(Xlag_diff), ncol = length(lags))
X <- c()

for (j in 1:ncol(Xlag_diff)) {
    for (i in seq_along(lags)) {
        x[, i] <- Hmisc::Lag(Xlag_diff[, j], lags[i])
    }
    X <- cbind(X, x)
}
