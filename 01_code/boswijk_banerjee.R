### local Parameters


Xlag <-cbind(q,r)

Y_dif <- diff(r)
W <- diff(Xlag)



if (lags > 0) {
    W <- lag(W, 1:lags)
    
}


