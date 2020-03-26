### Local Parameters
lags <- 2

Xlag <- cbind(depVar, indepVar)
Y_dif <- diff(depVar) # muss als numeric vorliegen
W <- diff(indepVar) # muss als numeric vorliegen

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

### Loop
res <- matrix(NA, nrow = nrow(Xlag) - lags - 1, ncol = ncol(Xlag))

for (i in 1:ncol(Xlag)) {
    loop_lm <- lm(Hmisc::Lag(Xlag[, i], shift = 1)[-1] ~ W)
    res[, i] <- as.numeric(loop_lm$residuals)
}

### Boswijk/Banerjee Test
BB_lm <- lm(Y_dif ~ W)
BB_res <- BB_lm$residuals

lm_res <- lm(BB_res ~ res[, 1] - res) # Interpretation des Minus?
betas <- coef(lm_res)
var <- vcov(lm_res)

stat[3] <- betas[1]/sqrt(var[1, 1]) # Welche Position des betas bei Stata?
stat[4] <- betas * solve(var) * betas
