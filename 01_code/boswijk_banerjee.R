### local Parameters


Xlag <- cbind(p,q)

Y_dif <- diff(r)
W <- diff(Xlag)


if (lags >= 1) {
    W <- lag(W)*diff(Xlag)
    
}


