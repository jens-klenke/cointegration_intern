# https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/PkgUtils


lags <- c(1,2,5)
df <- data.frame(
    x = c(11:20),
    y = c(21:30))

x <- matrix(NA, nrow = nrow(df), ncol =  length(lags))
x_names <- rep(NA, length(lags))
X <- c()
X_names <- NULL

for (j in 1:ncol(df)) {
    for (i in seq_along(lags)) {
    x[,i] <- Hmisc::Lag(df[,j], lags[i])
    x_names[i] <- paste0(names(df)[j],"_L_", lags[i] )
    }
    X <- cbind(X,x)
    X_names <- cbind(X_names, x_names)
    colnames(X) <- X_names
    }
