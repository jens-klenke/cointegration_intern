## packages
source(here::here('01_code/packages/packages.R'))

## sub functions 
source(here::here('09_simulation_and_approximation-cdf/sim_sub_functions.R'))


## parallelation 
num_cores <- detectCores() # need to change for the server
use_cores <- num_cores/2  # need to change for the server

registerDoParallel(use_cores) # cores need to be hard coded for the server

## Parameters
R2 <- seq(0, 0.95, 0.05) # long term correlation
N <- 1000 # length for every trajectory
c <- 0 # Parameter pesavanto model
K <- 11 #  max of Variables
rep <- 25000 # Number of Repetitions 25000
cases <- 3 # cases 
lambda <- (seq(1:N)/N) 

## matricies

Null_Distr_E <-  array(NA, dim = c(kmax, cases, rep) ) # matrix should be efficient 




## Loop

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
                if (cc == 1 && rr == 1){
                    BoswijkPValue <- 1 - rankindx(NullStatBoswijk,1)/rep+10^(-1000)
                }
                
                NullDistrJohansen <- sort(NullStatJohansen)
                if (cc == 1 && rr == 1){
                    JohansenPValue <- 1-rankindx(NullStatJohansen,1)/rep+10^(-1000)
                }
                
                NullDistrEngleGranger <- sort(NullStatEngleGranger)
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
                
                #4 tests
                FisherStat_B_ECR_J_E <- -2*(log(ErrCorrPValue)+log(BoswijkPValue)+log(EngleGrangerPValue)+log(JohansenPValue))
                
                
                if (cc == 1 && rr == 1){# Write Critical Values"
                    # 2 tests 
                    
                    # Engel - Johansen
                    NullDistrFisher_E_J <- sort(FisherStat_E_J)
                    
                    #4 tests
                    NullDistrFisher_B_ECR_J_E <- sort(FisherStat_B_ECR_J_E)
                    
                    
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
        

        # Engel -- Johansen
        Null_Distr_E_J[k,dets, ] <- NullDistrFisher_E_J
        
        #All Test
        Null_Distr_B_ECR_J_E[k,dets,]  <- NullDistrFisher_B_ECR_J_E
        
        
    }# dets 
}# k 


