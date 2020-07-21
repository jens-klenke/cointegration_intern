# Ohrnstein Uhlenbeck Process"
BU <- function(uu,d){
    tt <-dim(uu)[1] 
    rho <- (1+d/tt)
    v <- matrix(rep(0, length(uu)), nrow = tt)
    v[ ,1] <- uu[,1]
    for (t in 2:tt){
        v[t,] <- rho*v[t-1, ] + uu[t, ]
    }
    B <- v/sqrt(tt)
}

## rankindx
rankindx <- function(a,b){
    A <- dim(matrix(a))
    aux <-cbind(a, 1:A[1])
    aux <- aux[order(aux[,b]),]
    aux <- cbind(aux, 1:A[1])
    aux <- aux[order(aux[,A[2]+1]),]
    rankindx <- aux[ ,ncol(aux)]
    return(rankindx)
}