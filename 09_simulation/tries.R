T <- 10
k <- 3
R2run <- 0.05
lambda <- (seq(1:T)/T)

u <- matrix(c(0.0780,    1.9085,   -1.4158,    0.4147,
    1.3244,    0.1222,    0.0596,    0.3484,
    -0.2132,    1.0470,   -0.4113,    0.3493,
    -0.1345,   -0.2269,   -0.3680,   -0.7292,
    -1.1714,   -0.1625,   -1.3610,    0.3268,
    -1.3853,    0.6901,    0.7796,   -0.5149,
    0.3105,    0.5558,    0.4394,   -0.8964,
    -0.2495,   -1.1203,   -0.0896,   -1.2033,
    0.5037,   -1.5327,    1.0212,    1.0378,
    -0.8927,   -1.0979,   -0.8740,   -0.8459),
    ncol = k+1, byrow = TRUE)

#Zeile 198
W1 <- apply(matrix(u[, 1:k], ncol = k),2,cumsum)/sqrt(T)

#Zeile 199

u12 <- matrix(sqrt(R2run/(1-R2run))*(u[,1:k]%*% matrix(rep(1, k),ncol = 1))/sqrt(k) + u[ ,k+1], ncol = 1)

#Zeile 207
W1d <- W1 - (matrix(rep( apply(W1,2,mean), T), nrow = T, byrow = TRUE))

#zeile 209
W1 - (4-6*matrix(rep(t(lambda), k), ncol = k))*matrix(rep(apply(W1,2,mean),T), nrow = T, byrow = TRUE)- (12*matrix(rep(t(lambda), k), ncol = k)-6)*matrix(rep(apply(matrix(rep((lambda), k), ncol = k)*W1,2, mean), T), ncol = k, byrow = TRUE) 

#Zeile 210
J12dc <- J12 - (4-6*lambda)*rep(mean(J12), T) - (12*lambda-6) * rep(mean(lambda*J12), T)

#Zeile 218
WdcDW2 <- apply(Wdc[1:T-1,]* matrix(rep(u[2:T,k+1], k+1), ncol = k+1), 2, mean)

#Zeile 219
WdcWdci <- solve(1/T^2* t(Wdc)%*%Wdc)

#Zeile 220
W1dW1di <- solve(1/T * t(W1d[1:T-1,])%*% W1d[1:T-1,])

#Zeile 221
W1dJ12dc <- apply(W1d[1:T-1,] * matrix(rep(J12dc[1:T-1,],k),ncol = k), 2, mean) 

#Zeile 222
J12dc_sq <- mean(J12dc[1:T-1]^2)

#Zeile 223
J12DW2 <- mean(J12dc[1:T-1] * u[2:T,k+1])

## bosiwjk
#Zeile 225
BoswijkStat[j] <- c_run^2*J12dc_sq + 2*c_run*sqrt(T)*J12DW2 + WdcDW2%*%WdcWdci%*%WdcDW2

#Zeile 227
BoswijkPValue[j] <- 1- min(abs(BoswijkStat[j]- NullDistrBoswijk))/rep+10^(-10000)

#Zeile 232
Gc <- matrix(apply(Wdc*matrix(rep(J12dc, ncol(Wdc)), ncol = ncol(Wdc)),2 ,mean), ncol = 1)%*%t(matrix(c(rep(0,k), c_run)/sqrt(T), ncol = 1))



#Zeile 235
Wdc_dW_pr <- 1/T*t(u[2:T,])%*%Wdc[1:T-1,]

#Zeile 236
dW_Wdc_pr <- 1/T*t(Wdc[1:T-1,])%*%u[2:T,]

JohansenStat[j] <- max(eigen(Wdc_dW_pr%*%WdcWdci%*%dW_Wdc_pr+t(Gc)%*%WdcWdci%*%dW_Wdc_pr+ t(dW_Wdc_pr)%*%WdcWdci%*%Gc+t(Gc)%*%WdcWdci%*%Gc)$values)

#Zeile 245
etadc <-t(c( (-W1dW1di%*%(apply(W1d[1:T-1,] * matrix(rep(J12dc[1:T-1,], k),ncol = k), 2, mean))) , 1))

#Zeile 246
Adc <- 1/T*t(Wdc[1:T-1,])%*%Wdc[1:T-1,]

#Zeile 247
Dmat <- rbind(cbind(diag(k),  (sqrt(R2run/(1-R2run))*rep(1,k)/sqrt(k))), 
              c((sqrt(R2run/(1-R2run))*rep(1,k)/sqrt(k)), (1+R2run/(1-R2run))))

#Zeile 252
EngleGrangerStat[j] <- c_run*sqrt( t(etadc)%*%Adc%*%etadc)/sqrt(t(etadc)%*%Dmat%*%etadc)+
    (t(etadc)%*%Wdc_dWtilde%*%etadc)/(sqrt(t(etadc)%*%Dmat%*%etadc)*sqrt(t(etadc)%*%Adc%*%etadc))


#Zeile 262
zaehler <- sqrt(T)* (t(J12DW2)-W1dJ12dc%*%W1dW1di%*%apply(W1d[1:T-1,]*
                    matrix(rep(u[2:T,k+1],k), ncol = k), 2, mean ))

#Zeile 264
nenner <- (sqrt(as.complex( t(J12dc_sq) - W1dJ12dc%*%W1dW1di%*%W1dJ12dc)))

#Zeile 267
ErrCorrStat[j] <- c_run*sqrt(as.complex(t(J12dc_sq)-W1dJ12dc%*%W1dW1di%*%W1dJ12dc)) + zaehler/nenner

