## Ohrnstein Uhlenbeck Process
BU <- function(uu,d){
    tt <- dim(uu)[1] 
    rho <- (1+d/tt)
    v <- matrix(rep(0, length(uu)), nrow = tt)
    #v <- matrix(data = 0, nrow = tt, ncol = dim(uu)[2]) Eventuell minimal schneller
    v[, 1] <- uu[, 1]
    for (t in 2:tt){
        v[t, ] <- rho*v[t-1, ] + uu[t, ]
    }
    B <- v/sqrt(tt)
}

## rankindx
rankindx <- function(a, b){
    A <- dim(matrix(a)) 
    aux <- cbind(a, 1:A[1])
    aux <- aux[order(aux[, b]), ]
    aux <- cbind(aux, 1:A[1])
    aux <- aux[order(aux[, A[2]+1]), ]
    rankindx <- aux[, ncol(aux)]
    return(rankindx)
}

## term functions

sim_terms <- function(N, k, R2run, c_run, dets){
    u <- matrix(rnorm(N*(k+1)), nrow = N, ncol = k+1) # Draw random Shocks
    W1 <- apply(matrix(u[, 1:k], ncol = k), 2, cumsum)/sqrt(N)
    u12 <- matrix(sqrt(R2run/(1-R2run))*(u[, 1:k] %*% matrix(rep(1, k), ncol = 1))/sqrt(k) + u[, k+1], ncol = 1) 
    J12 <- BU(u12,c_run) # Ohrnstein Uhlenbeck Process
    
    # Corrections according to case (Pesavanto)
    
    if (dets == 1){ #No Constant, no trend"
        W1d <- W1
        J12dc <- J12
    } else if (dets == 2){ #Constant, no trend"
        W1d <- W1 - matrix(rep(apply(W1, 2, mean), N), nrow = N, byrow = TRUE)
        J12dc <- J12 - matrix(rep(mean(J12), N), nrow = N, byrow = TRUE)
       # J12dc <- J12 - matrix(rep(1, N)* apply(J12, 2, mean), nrow = N, byrow = TRUE)
    } else if (dets == 3){# Constant and Trend"
        W1d <- W1 - (4-6*matrix(rep(lambda, k), ncol = k))*matrix(rep(apply(W1, 2, mean), N), nrow = N, byrow = TRUE) - (12*matrix(rep(lambda, k), ncol = k)-6)*matrix(rep(apply(matrix(rep(lambda, k), ncol = k)*W1, 2, mean), N), ncol = k, byrow = TRUE) 
        # W1d <- W1 - (4-6*matrix(rep(t(lambda), k), ncol = k))*matrix(rep(apply(W1, 2, mean), N), nrow = N, byrow = TRUE) - (12*matrix(rep(t(lambda), k), ncol = k)-6)*matrix(rep(apply(matrix(rep((lambda), k), ncol = k)*W1,2, mean), N), ncol = k, byrow = TRUE) 
        J12dc <- J12 - (4-6*lambda)*rep(mean(J12), N) - (12*lambda-6) * rep(mean(lambda*J12), N) 
    }
    Wdc <- cbind(W1d, J12dc)
    
    # -----------------------------Common Terms------------"
    WdcDW2 <- apply(Wdc[1:N-1,]* matrix(rep(u[2:N,k+1], k+1), ncol = k+1), 2, mean)
    WdcWdci <- solve(1/N^2*t(Wdc)%*%Wdc)
    W1dW1di <- solve(1/N * t(W1d[1:N-1,])%*% W1d[1:N-1,])
    W1dJ12dc <- apply(W1d[1:N-1,] * matrix(rep(J12dc[1:N-1,],k),ncol = k), 2, mean) 
    J12dc_sq <- mean(J12dc[1:N-1]^2)
    J12DW2 <- mean(J12dc[1:N-1] * u[2:N,k+1])
    
    results  <- list(
        u = u,
        u12 = u12,
        W1d = W1d,
        Wdc = Wdc, 
        WdcDW2 = WdcDW2, 
        WdcWdci = WdcWdci, 
        W1dW1di = W1dW1di, 
        W1dJ12dc = W1dJ12dc,
        J12dc = J12dc,
        J12dc_sq = J12dc_sq, 
        J12DW2 = J12DW2) 
    
    return(results)    
}

# stat Boswijk
stat_Boswijk <- function(c_run, N, J12dc_sq, J12DW2, WdcDW2, WdcWdci){
    
    stat <- c_run^2*J12dc_sq + 2*c_run*sqrt(N)*J12DW2 + WdcDW2%*%WdcWdci%*%WdcDW2
    return(stat)

}

# stat Johansen

stat_Johansen <- function(Wdc, J12dc, u, N, k, c_run, WdcWdci){
    
    Gc <- matrix(apply(Wdc*matrix(rep(J12dc, ncol(Wdc)), ncol = ncol(Wdc)),2 ,mean), ncol = 1)%*%t(matrix(c(rep(0,k), c_run)/sqrt(N), ncol = 1))    
    Wdc_dW_pr <- 1/N*t(u[2:N,])%*%Wdc[1:N-1,]
    dW_Wdc_pr <- 1/N*t(Wdc[1:N-1,])%*%u[2:N,]
    
    stat <- max(eigen(Wdc_dW_pr%*%WdcWdci%*%dW_Wdc_pr+t(Gc)%*%WdcWdci%*%dW_Wdc_pr+ t(dW_Wdc_pr)%*%WdcWdci%*%Gc+t(Gc)%*%WdcWdci%*%Gc)$values)
    
    return(stat)
    
}

# stat Engel Granger

stat_Eng_Gr <- function(c_run, R2run, N, k, W1dW1di, W1d, J12dc, Wdc, u, u12){
    
    etadc <-matrix(c( (-W1dW1di%*%(apply(W1d[1:N-1,] * matrix(rep(J12dc[1:N-1,], k),ncol = k), 2, mean))) , 1), ncol = 1)
    Adc <- 1/N*t(Wdc[1:N-1,])%*%Wdc[1:N-1,]
    Dmat <- rbind(cbind(diag(k),  (sqrt(R2run/(1-R2run))*rep(1,k)/sqrt(k))), 
                  c((sqrt(R2run/(1-R2run))*rep(1,k)/sqrt(k)), (1+R2run/(1-R2run))))
    utilde <- cbind(u[,1:k], u12)
    Wdc_dWtilde <- 1/sqrt(N)*t(Wdc[1:N-1,])%*%utilde[2:N,]
    
    stat <- c_run*sqrt( t(etadc)%*%Adc%*%etadc)/sqrt(t(etadc)%*%Dmat%*%etadc)+
        (t(etadc)%*%Wdc_dWtilde%*%etadc)/(sqrt(t(etadc)%*%Dmat%*%etadc)*sqrt(t(etadc)%*%Adc%*%etadc))
    
    return(stat)
    
}

# ECR (Banerjee)

stat_Banerjee <- function(c_run, R2run, N, k, J12DW2, W1dJ12dc, W1dW1di, W1d, J12dc_sq, u){
    
    zaehler <- sqrt(N)* (t(J12DW2)-W1dJ12dc%*%W1dW1di%*%apply(W1d[1:N-1,]*
                                                                  matrix(rep(u[2:N,k+1],k), ncol = k), 2, mean ))               
    nenner <- (sqrt(as.complex( t(J12dc_sq) - W1dJ12dc%*%W1dW1di%*%W1dJ12dc)))
    
    stat <- c_run*sqrt(as.complex(t(J12dc_sq)-W1dJ12dc%*%W1dW1di%*%W1dJ12dc)) + zaehler/nenner
    
    stat <- Re(stat)
    
    return(stat)
    
}













