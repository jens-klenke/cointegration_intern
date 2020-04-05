#### Main part ####
#sub functions
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

###### Parameters
# savecv=0; # Save Critical values"
R2 <- seq(0, 0.95, 0.05)

T <- 10
#c <- c(-(seq(0,30,1)))
c = 0
kmax <- 11 #  max of Variables
rep <- 25000# Number of Repetitions 25000
cases <- 3 # cases 
lambda <- (seq(1:T)/T) # (1/T:1/T:1)





# Initalisierung 
# Engel
Null_Distr_E_J <-  array(NA, dim = c(kmax, cases, rep) )
cv_E_J <- matrix(NA, nrow = kmax, ncol = cases)
cv_E_J_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_E_J_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_E_J_3 <- matrix(NA, nrow = kmax, ncol = cases)

Null_Distr_E_B <-  array(NA, dim = c(kmax, cases, rep) )
cv_E_B <- matrix(NA, nrow = kmax, ncol = cases)
cv_E_B_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_E_B_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_E_B_3 <- matrix(NA, nrow = kmax, ncol = cases)

Null_Distr_E_ECR <-  array(NA, dim = c(kmax, cases, rep) )
cv_E_ECR <- matrix(NA, nrow = kmax, ncol = cases)
cv_E_ECR_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_E_ECR_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_E_ECR_3 <- matrix(NA, nrow = kmax, ncol = cases)

# Johansen
Null_Distr_J_E  <-  array(NA, dim = c(kmax, cases, rep) )
cv_J_E  <- matrix(NA, nrow = kmax, ncol = cases)
cv_J_E_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_J_E_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_J_E_3 <- matrix(NA, nrow = kmax, ncol = cases)

Null_Distr_J_B  <-  array(NA, dim = c(kmax, cases, rep) )
cv_J_B  <- matrix(NA, nrow = kmax, ncol = cases)
cv_J_B_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_J_B_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_J_B_3 <- matrix(NA, nrow = kmax, ncol = cases)

Null_Distr_J_ECR  <-  array(NA, dim = c(kmax, cases, rep) )
cv_J_ECR  <- matrix(NA, nrow = kmax, ncol = cases)
cv_J_ECR_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_J_ECR_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_J_ECR_3 <- matrix(NA, nrow = kmax, ncol = cases)

# Boswijk
Null_Distr_B_E  <-  array(NA, dim = c(kmax, cases, rep) )
cv_B_E  <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_E_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_E_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_E_3 <- matrix(NA, nrow = kmax, ncol = cases)

Null_Distr_B_J <-  array(NA, dim = c(kmax, cases, rep) )
cv_B_J  <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_J_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_J_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_J_3 <- matrix(NA, nrow = kmax, ncol = cases)


Null_Distr_B_ECR <-  array(NA, dim = c(kmax, cases, rep) )
cv_B_ECR  <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_3 <- matrix(NA, nrow = kmax, ncol = cases)

# Banerjee
Null_Distr_ECR_E <-  array(NA, dim = c(kmax, cases, rep) )
cv_ECR_E  <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_E_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_E_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_E_3 <- matrix(NA, nrow = kmax, ncol = cases)

Null_Distr_ECR_J <-  array(NA, dim = c(kmax, cases, rep) )
cv_ECR_J  <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_J_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_J_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_J_3 <- matrix(NA, nrow = kmax, ncol = cases)

Null_Distr_ECR_B <-  array(NA, dim = c(kmax, cases, rep) )
cv_ECR_B  <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_B_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_B_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_B_3 <- matrix(NA, nrow = kmax, ncol = cases)

# 3 Tests
Null_Distr_B_J_E <-  array(NA, dim = c(kmax, cases, rep) )
cv_ECR_B_J_E  <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_B_J_E_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_B_J_E_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_B_J_E_3 <- matrix(NA, nrow = kmax, ncol = cases)

Null_Distr_B_ECR_J <-  array(NA, dim = c(kmax, cases, rep) )
cv_B_ECR_J  <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_J_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_J_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_J_3 <- matrix(NA, nrow = kmax, ncol = cases)

Null_Distr_B_ECR_E <-  array(NA, dim = c(kmax, cases, rep) )
cv_B_ECR_E  <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_E_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_E_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_E_3 <- matrix(NA, nrow = kmax, ncol = cases)

Null_Distr_ECR_E_J <-  array(NA, dim = c(kmax, cases, rep) )
cv_ECR_E_J  <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_E_J_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_E_J_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_ECR_E_J_3 <- matrix(NA, nrow = kmax, ncol = cases)

#4 Test

Null_Distr_B_ECR_J_E <-  array(NA, dim = c(kmax, cases, rep) )
cv_B_ECR_J_E  <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_J_E_1 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_J_E_2 <- matrix(NA, nrow = kmax, ncol = cases)
cv_B_ECR_J_E_3 <- matrix(NA, nrow = kmax, ncol = cases)







### foreach -> k 

for (k in 1:kmax){# Number of Regressor Loop
    for (dets in 1:cases){# Number of cases Loop
        # initialization of Null Distribution of Test statistic
        NullStatBoswijk <- rep(NA, rep)
        NullStatJohansen <- rep(NA, rep)
        NullStatEngleGranger <- rep(NA, rep)
        NullStatErrCorr <- rep(NA, rep)
        NullDistrBoswijk <- NA
        NullDistrJohansen <- NA
        NullDistrEngleGranger <- NA
        NullDistrErrCorr <- NA

        for (rr in 1:length(R2)){# Loop over Pesavento R2"
            for (cc in 1:length(c)){#Loop over Local to unity parameter"
                # Initialize Test Statistics"
                BoswijkStat <- rep(NA, rep)
                JohansenStat <- rep(NA, rep)
                EngleGrangerStat <- rep(NA, rep)
                ErrCorrStat <- rep(NA, rep)
                BoswijkPValue <- rep(NA, rep)
                JohansenPValue <- rep(NA, rep)
                EngleGrangerPValue <- rep(NA, rep)
                ErrCorrPValue <- rep(NA, rep)
                
                # Set R2 and c corresponding to loop"
                R2run <- R2[rr]
                c_run <- c[cc]
                
                # Loop over repetitions"
                for (j in 1:rep){
                    u <- matrix(rnorm(T*(k+1)), nrow = T, ncol = k+1) # Draw random Shocks"
                    W1 <- apply(matrix(u[, 1:k], ncol = k),2,cumsum)/sqrt(T)
                    u12 <- matrix(sqrt(R2run/(1-R2run))*(u[,1:k]%*% matrix(rep(1, k),ncol = 1))/sqrt(k) + u[ ,k+1], ncol = 1)#Matrix multiplication
                    J12 <- BU(u12,c_run)#Ohrnstein Uhlenbeck Process"
                    # Corrections according to case"

                    if (dets==1){# No Constant, no trend"
                        W1d <- W1
                        J12dc <- J12
                    } else if (dets==2){#Constant, no trend"
                        W1d <- W1 - (matrix(rep( apply(W1,2,mean), T), nrow = T, byrow = TRUE)) 
                        J12dc <- J12 - matrix(rep(1, T)* apply(J12, 2, mean), nrow = T, byrow = TRUE) 
                    } else if (dets==3){# Constant and Trend"
                        W1d <- W1 - (4-6*matrix(rep(t(lambda), k), ncol = k))*matrix(rep(apply(W1, 2, mean), T), nrow = T, byrow = TRUE) - (12*matrix(rep(t(lambda), k), ncol = k)-6)*matrix(rep(apply(matrix(rep((lambda), k), ncol = k)*W1,2, mean), T), ncol = k, byrow = TRUE) 
                        J12dc <- J12 - (4-6*lambda)*rep(mean(J12), T) - (12*lambda-6) * rep(mean(lambda*J12), T) 
                    }
                    Wdc <- cbind(W1d, J12dc)

                    # -----------------------------Common Terms------------"
                    WdcDW2 <- apply(Wdc[1:T-1,]* matrix(rep(u[2:T,k+1], k+1), ncol = k+1), 2, mean)
                    WdcWdci <- solve(1/T^2* t(Wdc)%*%Wdc)
                    W1dW1di <- solve(1/T * t(W1d[1:T-1,])%*% W1d[1:T-1,])
                    W1dJ12dc <- apply(W1d[1:T-1,] * matrix(rep(J12dc[1:T-1,],k),ncol = k), 2, mean) 
                    J12dc_sq <- mean(J12dc[1:T-1]^2)
                    J12DW2 <- mean(J12dc[1:T-1] * u[2:T,k+1])
                    
                    # -------------------------------- Boswijk------------------------"
                    BoswijkStat[j] <- c_run^2*J12dc_sq + 2*c_run*sqrt(T)*J12DW2 + WdcDW2%*%WdcWdci%*%WdcDW2
                    if (cc == 1 && rr == 1){
                        NullStatBoswijk[j] <- BoswijkStat[j]
                    }
                    # -------------------------------- Johansen -------------------------------- */"

                    Gc <- matrix(apply(Wdc*matrix(rep(J12dc, ncol(Wdc)), ncol = ncol(Wdc)),2 ,mean), ncol = 1)%*%t(matrix(c(rep(0,k), c_run)/sqrt(T), ncol = 1))
                    Wdc_dW_pr <- 1/T*t(u[2:T,])%*%Wdc[1:T-1,]
                    dW_Wdc_pr <- 1/T*t(Wdc[1:T-1,])%*%u[2:T,]

                    JohansenStat[j] <- max(eigen(Wdc_dW_pr%*%WdcWdci%*%dW_Wdc_pr+t(Gc)%*%WdcWdci%*%dW_Wdc_pr+ t(dW_Wdc_pr)%*%WdcWdci%*%Gc+t(Gc)%*%WdcWdci%*%Gc)$values)
                    if (cc == 1 && rr == 1){
                        NullStatJohansen[j] <- JohansenStat[j]
                    }
                    
                    # -------------------------------- Engle-Granger -------------------------------- */"

                    etadc <-matrix(c( (-W1dW1di%*%(apply(W1d[1:T-1,] * matrix(rep(J12dc[1:T-1,], k),ncol = k), 2, mean))) , 1), ncol = 1)
                    Adc <- 1/T*t(Wdc[1:T-1,])%*%Wdc[1:T-1,]
                    Dmat <- rbind(cbind(diag(k),  (sqrt(R2run/(1-R2run))*rep(1,k)/sqrt(k))), 
                                  c((sqrt(R2run/(1-R2run))*rep(1,k)/sqrt(k)), (1+R2run/(1-R2run))))
                    utilde <- cbind(u[,1:k], u12)
                    Wdc_dWtilde <- 1/sqrt(T)*t(Wdc[1:T-1,])%*%utilde[2:T,]

                    EngleGrangerStat[j] <- c_run*sqrt( t(etadc)%*%Adc%*%etadc)/sqrt(t(etadc)%*%Dmat%*%etadc)+
                        (t(etadc)%*%Wdc_dWtilde%*%etadc)/(sqrt(t(etadc)%*%Dmat%*%etadc)*sqrt(t(etadc)%*%Adc%*%etadc))
                    
                    if (cc == 1 && rr == 1){
                        NullStatEngleGranger[j] <- EngleGrangerStat[j]
                    }

                    # -------------------------------- ECR (Banerjee) -------------------------------- */"
                    zaehler <- sqrt(T)* (t(J12DW2)-W1dJ12dc%*%W1dW1di%*%apply(W1d[1:T-1,]*
                               matrix(rep(u[2:T,k+1],k), ncol = k), 2, mean ))                    
                    nenner <- (sqrt(as.complex( t(J12dc_sq) - W1dJ12dc%*%W1dW1di%*%W1dJ12dc)))
                    ErrCorrStat[j] <- c_run*sqrt(as.complex(t(J12dc_sq)-W1dJ12dc%*%W1dW1di%*%W1dJ12dc)) + zaehler/nenner
                    
                    if (cc == 1){
                        NullStatErrCorr[j] <- Re(ErrCorrStat[j])
                    }

                }#rep loop end

                ## Write Null Distributions and Critical Values for underlying tests
                NullDistrBoswijk <- sort(NullStatBoswijk)
                CritvalBoswijk <- NullDistrBoswijk[rep*.95]
                CritvalBoswijk1 <- NullDistrBoswijk[rep*.99]
                CritvalBoswijk2 <- NullDistrBoswijk[rep*.95]
                CritvalBoswijk3 <- NullDistrBoswijk[rep*.90]
                if (cc == 1 && rr == 1){
                    BoswijkPValue <- 1 - rankindx(NullStatBoswijk,1)/rep+10^(-1000)
                }

                NullDistrJohansen <- sort(NullStatJohansen)
                CritvalJohansen <- NullDistrJohansen[rep*.95]
                CritvalJohansen1 <- NullDistrJohansen[rep*.99]
                CritvalJohansen2 <- NullDistrJohansen[rep*.95]
                CritvalJohansen3 <- NullDistrJohansen[rep*.90]
                if (cc == 1 && rr == 1){
                    JohansenPValue <- 1-rankindx(NullStatJohansen,1)/rep+10^(-1000)
                }

                NullDistrEngleGranger <- sort(NullStatEngleGranger)
                CritvalEngleGranger <- NullDistrEngleGranger[rep*0.05]
                CritvalEngleGranger1 <- NullDistrEngleGranger[rep*0.01]
                CritvalEngleGranger2 <- NullDistrEngleGranger[rep*0.05]
                CritvalEngleGranger3 <- NullDistrEngleGranger[rep*0.10]
                if (cc == 1 && rr == 1){
                    EngleGrangerPValue <- rankindx(NullStatEngleGranger,1)/rep+10^(-1000)
                }

                NullDistrErrCorr <- sort(NullStatErrCorr)
                CritvalErrCorr <- NullDistrErrCorr[rep*0.05]

                if (cc == 1){
                    ErrCorrPValue <- rankindx(NullStatErrCorr,1)/rep+10^(-1000)
                    if (rr==1){
                        CritValErrCorr1 <- NullDistrErrCorr[rep*0.01]
                        CritValErrCorr2 <- NullDistrErrCorr[rep*0.05]
                        CritValErrCorr3 <- NullDistrErrCorr[rep*0.10]
                    }
                }
         
## -------------------------------- Fisher Type Tests --------------------------------"
                # Define Statistics" # ErrCorrPValue
                # 2 tests 
                    # Engel
                FisherStat_E_J <- -2*(log(EngleGrangerPValue) + log(JohansenPValue))
                FisherStat_E_B <- -2*(log(EngleGrangerPValue) + log(BoswijkPValue))
                FisherStat_E_ECR <- -2*(log(EngleGrangerPValue) + log(ErrCorrPValue))
                    # JOhansen
                FisherStat_J_E <- -2*(log(JohansenPValue) + log(EngleGrangerPValue))
                FisherStat_J_B <- -2*(log(JohansenPValue) + log(BoswijkPValue))
                FisherStat_J_ECR <- -2*(log(JohansenPValue) + log(ErrCorrPValue))
                    # Boswijk
                FisherStat_B_E <- -2*(log(BoswijkPValue) + log(EngleGrangerPValue))
                FisherStat_B_J <- -2*(log(BoswijkPValue) + log(JohansenPValue))
                FisherStat_B_ECR <- -2*(log(BoswijkPValue) + log(ErrCorrPValue))
                    # Benerjee
                FisherStat_ECR_E <- -2*(log(ErrCorrPValue) + log(EngleGrangerPValue))
                FisherStat_ECR_J <- -2*(log(ErrCorrPValue) + log(JohansenPValue))
                FisherStat_ECR_B <- -2*(log(ErrCorrPValue) +log(BoswijkPValue))
                
                # 3 tests
                FisherStat_B_J_E <- -2*(log(EngleGrangerPValue)+log(JohansenPValue)+log(BoswijkPValue))
                FisherStat_B_ECR_J <- -2*(log(ErrCorrPValue)+log(BoswijkPValue)+log(JohansenPValue))
                FisherStat_B_ECR_E <- -2*(log(ErrCorrPValue)+log(BoswijkPValue)+log(EngleGrangerPValue)) # new Jens
                FisherStat_ECR_E_J <- -2*(log(ErrCorrPValue)+log(JohansenPValue)+log(EngleGrangerPValue)) # new Jens
                
                #4 tests
                FisherStat_B_ECR_J_E <- -2*(log(ErrCorrPValue)+log(BoswijkPValue)+log(EngleGrangerPValue)+log(JohansenPValue))
                
                
                FisherStatEJw=-4*((1/(1+exp(5*(R2run-.25))))*log(EngleGrangerPValue)+(1-1/(1+exp(5*(R2run-.25))))*log(JohansenPValue))
                
                if (cc == 1 && rr == 1){# Write Critical Values"
                    # 2 tests 
                    
                    # Engel 
                    NullDistrFisher_E_J <- sort(FisherStat_E_J)
                    CritvalFisher_E_J <- NullDistrFisher_E_J[rep*.95]
                    CritvalFisher_E_J_1 <- NullDistrFisher_E_J[rep*.99]
                    CritvalFisher_E_J_2 <- NullDistrFisher_E_J[rep*.95]
                    CritvalFisher_E_J_3 <- NullDistrFisher_E_J[rep*.90]
                    
                    NullDistrFisher_E_B <- sort(FisherStat_E_B)
                    CritvalFisher_E_B <- NullDistrFisher_E_B[rep*.95]
                    CritvalFisher_E_B_1 <- NullDistrFisher_E_B[rep*.99]
                    CritvalFisher_E_B_2 <- NullDistrFisher_E_B[rep*.95]
                    CritvalFisher_E_B_3 <- NullDistrFisher_E_B[rep*.90]
                    
                     
                    NullDistrFisher_E_ECR <- sort(FisherStat_E_ECR)
                    CritvalFisher_E_ECR <- NullDistrFisher_E_ECR[rep*.95]
                    CritvalFisher_E_ECR_1 <- NullDistrFisher_E_ECR[rep*.99]
                    CritvalFisher_E_ECR_2 <- NullDistrFisher_E_ECR[rep*.95]
                    CritvalFisher_E_ECR_3 <- NullDistrFisher_E_ECR[rep*.90]
                    
                    # Johansen
                    
                    NullDistrFisher_J_E <- sort(FisherStat_J_E)
                    CritvalFisher_J_E <- NullDistrFisher_J_E[rep*.95]
                    CritvalFisher_J_E_1 <- NullDistrFisher_J_E[rep*.99]
                    CritvalFisher_J_E_2 <- NullDistrFisher_J_E[rep*.95]
                    CritvalFisher_J_E_3 <- NullDistrFisher_J_E[rep*.90]
                    
                    NullDistrFisher_J_B <- sort(FisherStat_J_B)
                    CritvalFisher_J_B <- NullDistrFisher_J_E[rep*.95]
                    CritvalFisher_J_B_1 <- NullDistrFisher_J_B[rep*.99]
                    CritvalFisher_J_B_2 <- NullDistrFisher_J_B[rep*.95]
                    CritvalFisher_J_B_3 <- NullDistrFisher_J_B[rep*.90]
                    
                    NullDistrFisher_J_ECR <- sort(FisherStat_J_ECR) 
                    CritvalFisher_J_ECR <- NullDistrFisher_J_ECR[rep*.95]
                    CritvalFisher_J_ECR_1 <- NullDistrFisher_J_ECR[rep*.99]
                    CritvalFisher_J_ECR_2 <- NullDistrFisher_J_ECR[rep*.95]
                    CritvalFisher_J_ECR_3 <- NullDistrFisher_J_ECR[rep*.90]
                    
                    # Boswijk

                    NullDistrFisher_B_E <- sort(FisherStat_B_E)
                    CritvalFisher_B_E <- NullDistrFisher_B_E[rep*.95]
                    CritvalFisher_B_E_1 <- NullDistrFisher_B_E[rep*.99]
                    CritvalFisher_B_E_2 <- NullDistrFisher_B_E[rep*.95]
                    CritvalFisher_B_E_3 <- NullDistrFisher_B_E[rep*.90]
                    
                    NullDistrFisher_B_J <- sort(FisherStat_B_J)
                    CritvalFisher_B_J <- NullDistrFisher_B_J[rep*.95]
                    CritvalFisher_B_J_1 <- NullDistrFisher_B_J[rep*.99]
                    CritvalFisher_B_J_2 <- NullDistrFisher_B_J[rep*.95]
                    CritvalFisher_B_J_3 <- NullDistrFisher_B_J[rep*.90]

                    NullDistrFisher_B_ECR <- sort(FisherStat_B_ECR)
                    CritvalFisher_B_ECR <- NullDistrFisher_B_J[rep*.95]
                    CritvalFisher_B_ECR_1 <- NullDistrFisher_B_J[rep*.99]
                    CritvalFisher_B_ECR_2 <- NullDistrFisher_B_J[rep*.95]
                    CritvalFisher_B_ECR_3 <- NullDistrFisher_B_J[rep*.90]
                    
                    # Banerjee
                    
                    NullDistrFisher_ECR_E <- sort(FisherStat_ECR_E)
                    CritvalFisher_ECR_E <- NullDistrFisher_ECR_E[rep*.95]
                    CritvalFisher_ECR_E_1 <- NullDistrFisher_ECR_E[rep*.99]
                    CritvalFisher_ECR_E_2 <- NullDistrFisher_ECR_E[rep*.95]
                    CritvalFisher_ECR_E_3 <- NullDistrFisher_ECR_E[rep*.90]
                    
                    NullDistrFisher_ECR_J <- sort(FisherStat_ECR_J)
                    CritvalFisher_ECR_J <- NullDistrFisher_ECR_J[rep*.95]
                    CritvalFisher_ECR_J_1 <- NullDistrFisher_ECR_J[rep*.99]
                    CritvalFisher_ECR_J_2 <- NullDistrFisher_ECR_J[rep*.95]
                    CritvalFisher_ECR_J_3 <- NullDistrFisher_ECR_J[rep*.90]
                    
                    NullDistrFisher_ECR_B <- sort(FisherStat_ECR_B)
                    CritvalFisher_ECR_B <- NullDistrFisher_ECR_B[rep*.95]
                    CritvalFisher_ECR_B_1 <- NullDistrFisher_ECR_B[rep*.99]
                    CritvalFisher_ECR_B_2 <- NullDistrFisher_ECR_B[rep*.95]
                    CritvalFisher_ECR_B_3 <- NullDistrFisher_ECR_B[rep*.90]
                    
                    
                    
                    # 3 Tests
                    NullDistrFisher_B_J_E <- sort(FisherStat_B_J_E)
                    CritvalFisher_B_J_E <- NullDistrFisher_B_J_E[rep*.95]
                    CritvalFisher_B_J_E_1 <- NullDistrFisher_B_J_E[rep*.99]
                    CritvalFisher_B_J_E_2 <- NullDistrFisher_B_J_E[rep*.95]
                    CritvalFisher_B_J_E_3 <- NullDistrFisher_B_J_E[rep*.90]

                    NullDistrFisher_B_ECR_J <- sort(FisherStat_B_ECR_J)
                    CritvalFisher_B_ECR_J <- NullDistrFisher_B_ECR_J[rep*.95]
                    CritvalFisher_B_ECR_J_1 <- NullDistrFisher_B_ECR_J[rep*.99]
                    CritvalFisher_B_ECR_J_2 <- NullDistrFisher_B_ECR_J[rep*.95]
                    CritvalFisher_B_ECR_J_3 <- NullDistrFisher_B_ECR_J[rep*.90]
                    
                    NullDistrFisher_B_ECR_E <- sort(FisherStat_B_ECR_E)
                    CritvalFisher_B_ECR_E <- NullDistrFisher_B_ECR_E[rep*.95]
                    CritvalFisher_B_ECR_E_1 <- NullDistrFisher_B_ECR_E[rep*.99]
                    CritvalFisher_B_ECR_E_2 <- NullDistrFisher_B_ECR_E[rep*.95]
                    CritvalFisher_B_ECR_E_3 <- NullDistrFisher_B_ECR_E[rep*.90]
                    
                    NullDistrFisher_ECR_E_J <- sort(FisherStat_ECR_E_J)
                    CritvalFisher_ECR_E_J <- NullDistrFisher_ECR_E_J[rep*.95]
                    CritvalFisher_ECR_E_J_1 <- NullDistrFisher_ECR_E_J[rep*.99]
                    CritvalFisher_ECR_E_J_2 <- NullDistrFisher_ECR_E_J[rep*.95]
                    CritvalFisher_ECR_E_J_3 <- NullDistrFisher_ECR_E_J[rep*.90]
                    

                    #4 tests
                    NullDistrFisher_B_ECR_J_E <- sort(FisherStat_B_ECR_J_E)
                    CritvalFisher_B_ECR_J_E <- NullDistrFisher_B_ECR_J_E[rep*.95]
                    CritvalFisher_B_ECR_J_E_1 <- NullDistrFisher_B_ECR_J_E[rep*.99]
                    CritvalFisher_B_ECR_J_E_2 <- NullDistrFisher_B_ECR_J_E[rep*.95]
                    CritvalFisher_B_ECR_J_E_3 <- NullDistrFisher_B_ECR_J_E[rep*.90]
                    
                    
                    
               }

                print('------------------')
                print(paste('Simultet: R2 = ', R2run, 'which is the', rr ,'of' , length(R2)))
                print('------------------')
                print('')
                print(paste('Simultet: Case = ', dets, 'which is the', dets ,'of' , max(cases)))
                print('')
                print('------------------')
                print('')
                

            }# R2"
        }# c"
        
        
###### storing 
        
        #Füllen
        # Engel
        Null_Distr_E_J[k,dets, ] <- NullDistrFisher_E_J
        cv_E_J[k,dets] <- CritvalFisher_E_J
        cv_E_J_1[k,dets] <- CritvalFisher_E_J_1
        cv_E_J_2[k,dets] <- CritvalFisher_E_J_2
        cv_E_J_3[k,dets] <- CritvalFisher_E_J_3
        
        Null_Distr_E_B[k,dets,] <- NullDistrFisher_E_B
        cv_E_B[k,dets] <- CritvalFisher_E_B
        cv_E_B_1[k,dets] <- CritvalFisher_E_B_1
        cv_E_B_2[k,dets] <- CritvalFisher_E_B_2
        cv_E_B_3[k,dets] <- CritvalFisher_E_B_3
        
        Null_Distr_E_ECR[k,dets,] <- NullDistrFisher_E_ECR 
        cv_E_ECR[k,dets] <- CritvalFisher_E_ECR
        cv_E_ECR_1[k,dets] <- CritvalFisher_E_ECR_1
        cv_E_ECR_2[k,dets] <- CritvalFisher_E_ECR_2 
        cv_E_ECR_3[k,dets] <- CritvalFisher_E_ECR_3
        
        # Johansen
        Null_Distr_J_E[k,dets,] <- NullDistrFisher_J_E
        cv_J_E[k,dets]  <- CritvalFisher_J_E
        cv_J_E_1[k,dets] <- CritvalFisher_J_E_1
        cv_J_E_2[k,dets] <- CritvalFisher_J_E_2
        cv_J_E_3[k,dets] <- CritvalFisher_J_E_3
        
        Null_Distr_J_B[k,dets,]  <- NullDistrFisher_J_B
        cv_J_B[k,dets]  <- CritvalFisher_J_B
        cv_J_B_1[k,dets] <- CritvalFisher_J_B_1
        cv_J_B_2[k,dets] <- CritvalFisher_J_B_2
        cv_J_B_3[k,dets] <- CritvalFisher_J_B_3
        
        Null_Distr_J_ECR[k,dets,]  <-  NullDistrFisher_J_ECR
        cv_J_ECR[k,dets]  <- CritvalFisher_J_ECR
        cv_J_ECR_1[k,dets] <- CritvalFisher_J_ECR_1
        cv_J_ECR_2[k,dets] <- CritvalFisher_J_ECR_2
        cv_J_ECR_3[k,dets] <- CritvalFisher_J_ECR_3
        
        # Boswijk
        Null_Distr_B_E[k,dets,]  <- NullDistrFisher_B_E 
        cv_B_E[k,dets]  <- CritvalFisher_B_E
        cv_B_E_1[k,dets] <- CritvalFisher_B_E_1
        cv_B_E_2[k,dets] <- CritvalFisher_B_E_2
        cv_B_E_3[k,dets] <- CritvalFisher_B_E_3
        
        Null_Distr_B_J[k,dets,] <- NullDistrFisher_B_J
        cv_B_J[k,dets]  <- CritvalFisher_B_J
        cv_B_J_1[k,dets] <- CritvalFisher_B_J_1
        cv_B_J_2[k,dets] <- CritvalFisher_B_J_2
        cv_B_J_3[k,dets] <- CritvalFisher_B_J_3
        
        Null_Distr_B_ECR[k,dets,] <- NullDistrFisher_B_ECR
        cv_B_ECR[k,dets]  <- CritvalFisher_B_ECR
        cv_B_ECR_1[k,dets] <- CritvalFisher_B_ECR_1
        cv_B_ECR_2[k,dets] <- CritvalFisher_B_ECR_2
        cv_B_ECR_3[k,dets] <- CritvalFisher_B_ECR_3
        
        #Banerjee
        Null_Distr_ECR_E[k,dets,] <- NullDistrFisher_ECR_E
        cv_ECR_E[k,dets]  <- CritvalFisher_ECR_E
        cv_ECR_E_1[k,dets] <- CritvalFisher_ECR_E_1
        cv_ECR_E_2[k,dets] <- CritvalFisher_ECR_E_2
        cv_ECR_E_3[k,dets] <- CritvalFisher_ECR_E_3
        
        Null_Distr_ECR_J[k,dets,] <- NullDistrFisher_ECR_J
        cv_ECR_J[k,dets]  <- CritvalFisher_ECR_J
        cv_ECR_J_1[k,dets] <- CritvalFisher_ECR_J_1
        cv_ECR_J_2[k,dets] <- CritvalFisher_ECR_J_2
        cv_ECR_J_3[k,dets] <- CritvalFisher_ECR_J_3
        
        Null_Distr_ECR_B[k,dets,] <- NullDistrFisher_ECR_B
        cv_ECR_B[k,dets]  <- CritvalFisher_ECR_B
        cv_ECR_B_1[k,dets] <- CritvalFisher_ECR_B_1
        cv_ECR_B_2[k,dets] <- CritvalFisher_ECR_B_2
        cv_ECR_B_3[k,dets] <- CritvalFisher_ECR_B_3
        
        ### 3. Tests 
        
        Null_Distr_B_J_E[k,dets,] <- NullDistrFisher_B_J_E
        cv_ECR_B_J_E[k,dets]  <- CritvalFisher_B_J_E
        cv_ECR_B_J_E_1[k,dets] <- CritvalFisher_B_J_E_1
        cv_ECR_B_J_E_2[k,dets] <- CritvalFisher_B_J_E_2
        cv_ECR_B_J_E_3[k,dets] <- CritvalFisher_B_J_E_3
        
        
        Null_Distr_B_ECR_J[k,dets,] <- NullDistrFisher_B_ECR_J 
        cv_B_ECR_J[k,dets]  <- CritvalFisher_B_ECR_J
        cv_B_ECR_J_1[k,dets] <- CritvalFisher_B_ECR_J_1
        cv_B_ECR_J_2[k,dets] <- CritvalFisher_B_ECR_J_2
        cv_B_ECR_J_3[k,dets] <- CritvalFisher_B_ECR_J_3
        
        Null_Distr_B_ECR_E[k,dets,]  <- NullDistrFisher_B_ECR_E 
        cv_B_ECR_E[k,dets]   <- CritvalFisher_B_ECR_E
        cv_B_ECR_E_1[k,dets] <- CritvalFisher_B_ECR_E_1
        cv_B_ECR_E_2[k,dets] <- CritvalFisher_B_ECR_E_2
        cv_B_ECR_E_3[k,dets] <- CritvalFisher_B_ECR_E_3
        
        Null_Distr_ECR_E_J[k,dets,] <-  NullDistrFisher_ECR_E_J
        cv_ECR_E_J[k,dets]  <- CritvalFisher_ECR_E_J
        cv_ECR_E_J_1[k,dets] <- CritvalFisher_ECR_E_J_1
        cv_ECR_E_J_2[k,dets] <- CritvalFisher_ECR_E_J_2
        cv_ECR_E_J_3[k,dets] <- CritvalFisher_ECR_E_J_3
        
        #4 Test
        
        Null_Distr_B_ECR_J_E[k,dets,]  <- NullDistrFisher_B_ECR_J_E
        cv_B_ECR_J_E[k,dets]   <- CritvalFisher_B_ECR_J_E
        cv_B_ECR_J_E_1[k,dets] <- CritvalFisher_B_ECR_J_E_1
        cv_B_ECR_J_E_2[k,dets] <- CritvalFisher_B_ECR_J_E_2
        cv_B_ECR_J_E_3[k,dets] <- CritvalFisher_B_ECR_J_E_3

        
    }# dets */"
}# k */"



