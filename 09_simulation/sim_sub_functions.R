# Ohrnstein Uhlenbeck Process"
BU <- function(u,d){
tt <- nrow(u) ## glaube unnötig
rho <- (1+d/tt)
v <- matrix(rep(0, length(u)), nrow = tt)
v[ ,1] <- u[,1]
for (t in 2:tt){
    v[t,] <- rho*v[t-1, ] + u[t, ]
    }
B <- v/sqrt(tt)
}

###rankindx
rankindx <- function(a,b){
A <- dim(a)
aux <-cbind(a, 1:A[1])
aux <- aux[order(aux[,b]),]
aux <- cbind(aux, 1:A[1])
aux <- aux[order(aux[,A[2]+1]),]
rankindx <- aux[ ,ncol(aux)]
return(rankindx)
}

