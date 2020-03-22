# https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/PkgUtils


lags <- 2
df <- matrix(seq(1:20), ncol = 2)

x <- matrix(NA, nrow = nrow(df), ncol =  lags) 
X <- c()

for (j in 1:ncol(df)) {
    for (i in 1:lags) {
    x[,i] <- Hmisc::Lag(df[,j], i)    
    }
    X <- cbind(X,x)    
}
 
X
