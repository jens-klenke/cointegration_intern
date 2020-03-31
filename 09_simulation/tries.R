T <- 10
k <- 1

u <- matrix(c(    0.5201,   -0.2938,
                  -0.0200,   -0.8479,
                  -0.0348,   -1.1201,
                  -0.7982,    2.5260,
                  1.0187,    1.6555,
                  -0.1332,    0.3075,
                  -0.7145,   -1.2571,
                  1.3514,   -0.8655,
                  -0.2248,   -0.1765,
                  -0.5890,    0.7914),
    ncol = k+1, byrow = TRUE)

#Zeile 198
W1 <- apply(matrix(u[, 1:k], ncol = k),2,cumsum)/sqrt(T)

#Zeile 206
W1d <- W1 - matrix(rep( apply(W1,2,mean), T), nrow = T, byrow = TRUE)

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
Gc <- matrix(kronecker(apply(Wdc*matrix(rep(J12dc, nrow(Wdc)), ncol = length(J12dc)),2 ,mean),c(rep(0,k), c_run)/sqrt(T)), ncol = k+1, byrow = TRUE)# the square root is to make it fit with other extra power, see notes */"

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

EngleGrangerStat[j] <- c_run*sqrt( etadc%*%Adc%*%t(etadc))/sqrt(etadc%*%Dmat%*%t(etadc))+
    (etadc%*%Wdc_dWtilde%*%t(etadc))/(sqrt(etadc%*%Dmat%*%t(etadc))*sqrt(etadc%*%Adc%*%t(etadc)))

#Zeile 262
zaehler <- sqrt(T)* (t(J12DW2)-W1dJ12dc%*%W1dW1di%*%apply(W1d[1:T-1,]*
                    matrix(rep(u[2:T,k+1],k), ncol = k), 2, mean ))

#Zeile 264
nenner <- (sqrt(as.complex( t(J12dc_sq) - W1dJ12dc%*%W1dW1di%*%W1dJ12dc)))

#Zeile 267
ErrCorrStat[j] <- c_run*sqrt(as.complex(t(J12dc_sq)-W1dJ12dc%*%W1dW1di%*%W1dJ12dc)) + zaehler/nenner

