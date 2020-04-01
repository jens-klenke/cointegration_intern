#### Parameters
#load(here::here('09_simulation/functions.R'))
#asy <- function(savecv){
# This function calculates/simulates the asymptotic distributions of the
# various co-integration tests. It returns a
# scalar equal to 1 when done.
# If SAVECV is equal to 1, it saves the Array of Critical values CV
# in critical_values.mat and the distributions under the Null
# NULLDISTR in NullDistributions.mat
# In addition it provides asymptotic power results, results saved in
# workspace.mat
#--------------------------------------------------------------------------
# The Array NullDistr is organized as follows (neccessary for asytest.m)
# NullDistr cols  1: EngleGranger,      CV dim 3: nvars-1
#                 2: Johansen,          CV dim 4: deterministics (trend+1)
#                 3: ErrCorr,
#                 4: Boswijk
#--------------------------------------------------------------------------
# Critical values in critval.mat, Array CV are organised as follows:
#  CV: rows:  1:  Boswijk           CV cols : nvars-1
#             2:  ErrCorr           CV dim 3: deterministics (trend+1)
#             3:  Johansen          CV dim 4: LvL 1# 5# 10#
#             4:  EngleGranger
#             5:  FisherBECR;
#             6:  FisherBJ
#             7:  FisherBE
#             8:  FisherEJ
#             9:  FisherBJE
#             10: FisherBECRJ
#             11: FisherBECRJE"
#             12: InvNormBJE
#             13: MinBECR
#             14: MinEJ
#--------------------------------------------------------------------------"

#--------------------------------------------------------------------------"

#asy(1) # function call


## UR-Test Correction factors (Obtained in different GAUSS Code)
# UR-EJ"

 CritValCorrJohansenMultEJ <-  matrix( c(1.1, 1.077, 1.065,
                                        1.08, 1.076, 1.068,
                                        1.074, 1.063, 1.064,
                                        1.066, 1.059, 1.056,
                                        1.061, 1.055, 1.053,
                                        1.052, 1.051, 1.052,
                                        1.049, 1.047, 1.054,
                                        1.045, 1.045, 1.043,
                                        1.045, 1.042, 1.043,
                                        1.043, 1.043, 1.038,
                                        1.04, 1.039, 1.037),
                                      ncol = 3, byrow = TRUE)

CritValCorrEngleGrangerMultEJ <-  matrix( c( 1.065, 1.05, 1.043,
                                             1.058, 1.052, 1.044,
                                             1.055, 1.049, 1.046,
                                             1.051,	1.045, 1.042,
                                             1.048,	1.045, 1.041,
                                             1.046, 1.044, 1.04,
                                             1.045, 1.042, 1.035,
                                             1.042, 1.041, 1.039,
                                             1.04,	1.038, 1.039,
                                             1.039,	1.035, 1.037,
                                             1.038,	1.037, 1.035),
                                          ncol = 3, byrow = TRUE)

# UR-BJ

CritValCorrJohansenMultBJ <- matrix( c( 1.101, 1.083, 1.07,
                                        1.084, 1.082, 1.075,
                                        1.075, 1.067, 1.068,
                                        1.071, 1.063, 1.061,
                                        1.063, 1.058, 1.055,
                                        1.056, 1.052, 1.054,
                                        1.05, 1.053, 1.049,
                                        1.047, 1.048, 1.045,
                                        1.044, 1.042, 1.046,
                                        1.044, 1.161, 1.039,
                                        1.043, 1.039, 1.039),
                                     ncol = 3 , byrow = TRUE)

CritValCorrBoswijkMultBJ <-  matrix( c( 1.128, 1.104, 1.093,
                                        1.131, 1.11, 1.095,
                                        1.122, 1.104, 1.096,
                                        1.107, 1.099, 1.09,
                                        1.103, 1.094, 1.088,
                                        1.096, 1.091, 1.085,
                                        1.092, 1.082, 1.082,
                                        1.089, 1.08, 1.081,
                                        1.085, 1.081, 1.078,
                                        1.079, 1.008, 1.075,
                                        1.072, 1.076, 1.071),
                                     ncol = 3, byrow = TRUE)

# UR-BERC
CritValCorrErrCorrMultBERC <- matrix( c( 1.049, 1.022, 1.018,
                                         1.046,	1.028, 1.023,
                                         1.046, 1.033, 1.023,
                                         1.042, 1.033, 1.028,
                                         1.04, 1.032, 1.029,
                                         1.041, 1.034, 1.028,
                                         1.039, 1.035, 1.029,
                                         1.036, 1.032, 1.028,
                                         1.034,	1.032, 1.028,
                                         1.034, 1.031, 1.03,
                                         1.035, 1.032, 1.028),
                                      ncol = 3, byrow = TRUE)

CritValCorrBoswijkMultBERC <- matrix( c( 1.077, 1.042, 1.032,
                                         1.075, 1.052, 1.038,
                                         1.07, 1.053, 1.038,
                                         1.057, 1.053, 1.043,
                                         1.058, 1.049, 1.043,
                                         1.06,	1.051, 1.044,
                                         1.056,	1.055, 1.045,
                                         1.05,	1.044, 1.044,
                                         1.049,	1.047, 1.044,
                                         1.046,	1.041, 1.043,
                                         1.047,	1.045, 1.041),
                                      ncol = 3, byrow = TRUE)

#### Main part ####
#load sim_sub_function by hand!!
 
# savecv=0; # Save Critical values"
R2 <- seq(0, 0.95, 0.05)
# R2=0
T <- 100
#c <- c(-(seq(0,30,1)))
c = 0
kmax <- 3 #  max of Variables
rep <- 25 # Number of Repetitions 25000
cases <- 3 # cases 
lambda <- (seq(1:T)/T) # (1/T:1/T:1)

# Initialize Power Curves:"
# (1) Underlying Tests
BoswijkLocalAsyPower <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
JohansenLocalAsyPower <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
EngleGrangerLocalAsyPower <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
ErrCorrLocalAsyPower <-array(NA, dim = c(kmax, cases, length(c), length(R2)))

#(2) UR-Tests
UnionRejectionLocalAsyPowerBJAsymmMult <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
UnionRejectionLocalAsyPowerBERCAsymmMult <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
UnionRejectionLocalAsyPowerEJAsymmMult <- array(NA, dim = c(kmax, cases, length(c), length(R2)))

#(3) Fisher type Tests
FisherLocalAsyPowerEJ <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
FisherLocalAsyPowerBJ <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
FisherLocalAsyPowerBE <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
FisherLocalAsyPowerBECR <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
FisherLocalAsyPowerBJE <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
FisherLocalAsyPowerEJw <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
FisherLocalAsyPowerBECRJE <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
FisherLocalAsyPowerBECRJ <- array(NA, dim = c(kmax, cases, length(c), length(R2)))

#(4) Inverse Normal type Test
#InvNormLocalAsyPowerBJE <- array(NA, dim = c(kmax, cases, length(c), length(R2)))

#(5) Minimum tests"
MinLocalAsyPowerBERC <- array(NA, dim = c(kmax, cases, length(c), length(R2)))
MinLocalAsyPowerEJ <- array(NA, dim = c(kmax, cases, length(c), length(R2)))

for (k in 1:kmax){# Number of Regressor Loop"
    for (dets in 1:cases){# Number of cases Loop"
        # initialization of Null Distribution of Test statistic"
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
                    u <- matrix( rnorm(T*(k+1)), nrow = T, ncol = k+1) # Draw random Shocks"
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
                    } else {
                        BoswijkPValue[j] <- 1- min(abs(BoswijkStat[j]- NullDistrBoswijk))/rep+10^(-10000)
                    }

                    # -------------------------------- Johansen -------------------------------- */"

                    Gc <- matrix(apply(Wdc*matrix(rep(J12dc, ncol(Wdc)), ncol = ncol(Wdc)),2 ,mean), ncol = 1)%*%t(matrix(c(rep(0,k), c_run)/sqrt(T), ncol = 1))
                    Wdc_dW_pr <- 1/T*t(u[2:T,])%*%Wdc[1:T-1,]
                    dW_Wdc_pr <- 1/T*t(Wdc[1:T-1,])%*%u[2:T,]

                    JohansenStat[j] <- max(eigen(Wdc_dW_pr%*%WdcWdci%*%dW_Wdc_pr+t(Gc)%*%WdcWdci%*%dW_Wdc_pr+ t(dW_Wdc_pr)%*%WdcWdci%*%Gc+t(Gc)%*%WdcWdci%*%Gc)$values)
                    if (cc == 1 && rr == 1){
                        NullStatJohansen[j] <- JohansenStat[j]
                    } else {
                        JohansenPValue [j] <- 1-min(abs(JohansenStat[j]-NullDistrJohansen))/rep+10.^(-1000)
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
                    } else {
                        EngleGrangerPValue[j] <- min(abs(EngleGrangerStat[j]-NullDistrEngleGranger))/rep
                    }

                    # -------------------------------- ECR (Banerjee) -------------------------------- */"
                    zaehler <- sqrt(T)* (t(J12DW2)-W1dJ12dc%*%W1dW1di%*%apply(W1d[1:T-1,]*
                               matrix(rep(u[2:T,k+1],k), ncol = k), 2, mean ))                    
                    nenner <- (sqrt(as.complex( t(J12dc_sq) - W1dJ12dc%*%W1dW1di%*%W1dJ12dc)))
                    ErrCorrStat[j] <- c_run*sqrt(as.complex(t(J12dc_sq)-W1dJ12dc%*%W1dW1di%*%W1dJ12dc)) + zaehler/nenner
                    
                    if (cc == 1){# here we condition on each R2 because the r2=0 c.v. turned out not to be good for R2 very large >0.75*/"
                        NullStatErrCorr[j] <- Re(ErrCorrStat[j])# a 0 plus imaginary occurred once... */"
                    } else {
                        ErrCorrPValue[j] <- min(abs(Re(ErrCorrStat[j])-NullDistrErrCorr))/rep
                    }

                }#rep loop end"

                ## Write Null Distributions and Critical Values for underlying tests"
                NullDistrBoswijk <- sort(NullStatBoswijk)
                CritvalBoswijk <- NullDistrBoswijk[rep*.95]
                CritvalBoswijk1 <- NullDistrBoswijk[rep*.99]
                CritvalBoswijk2 <- NullDistrBoswijk[rep*.95]
                CritvalBoswijk3 <- NullDistrBoswijk[rep*.90]
                if (cc == 1 && rr == 1){
                    BoswijkPValue <- 1- rankindx(NullStatBoswijk,1)/rep+10.^(-1000)
                }
                BoswijkLocalAsyPower[k,dets,cc,rr] <- mean((BoswijkStat > CritvalBoswijk))

                NullDistrJohansen <- sort(NullStatJohansen)
                CritvalJohansen <- NullDistrJohansen[rep*.95]
                CritvalJohansen1 <- NullDistrJohansen[rep*.99]
                CritvalJohansen2 <- NullDistrJohansen[rep*.95]
                CritvalJohansen3 <- NullDistrJohansen[rep*.90]
                if (cc == 1 && rr == 1){
                    JohansenPValue <- 1-rankindx(NullStatJohansen,1)/rep+10.^(-1000)
                }
                JohansenLocalAsyPower[k,dets,cc,rr] <- mean((JohansenStat > CritvalJohansen))

                NullDistrEngleGranger <- sort(NullStatEngleGranger)
                CritvalEngleGranger <- NullDistrEngleGranger[rep*0.05]
                CritvalEngleGranger1 <- NullDistrEngleGranger[rep*0.01]
                CritvalEngleGranger2 <- NullDistrEngleGranger[rep*0.05]
                CritvalEngleGranger3 <- NullDistrEngleGranger[rep*0.10]
                if (cc == 1 && rr == 1){
                    EngleGrangerPValue <- rankindx(NullStatEngleGranger,1)/rep+10.^(-1000)
                }
                EngleGrangerLocalAsyPower[k,dets,cc,rr] = mean((EngleGrangerStat <= CritvalEngleGranger))

                NullDistrErrCorr <- sort(NullStatErrCorr)
                CritvalErrCorr <- NullDistrErrCorr[rep*0.05]

                if (cc == 1){
                    ErrCorrPValue <- rankindx(NullStatErrCorr,1)/rep+10.^(-1000)
                    if (rr==1){
                        CritValErrCorr1 <- NullDistrErrCorr[rep*0.01]
                        CritValErrCorr2 <- NullDistrErrCorr[rep*0.05]
                        CritValErrCorr3 <- NullDistrErrCorr[rep*0.10]
                    }
                }
                ErrCorrLocalAsyPower[k,dets,cc,rr] = mean((Re(ErrCorrStat) <= Re(CritvalErrCorr)))

                # -------------------------------- Union of Rejections -------------------------------- */"
                UnionRejectionStatEJAsymmMult <- (EngleGrangerStat*(EngleGrangerStat < CritvalEngleGranger*CritValCorrEngleGrangerMultEJ[k,dets]) < CritvalEngleGranger*CritValCorrEngleGrangerMultEJ[k,dets]) + (JohansenStat*(EngleGrangerStat > CritvalEngleGranger*CritValCorrEngleGrangerMultEJ[k,dets]) > CritvalJohansen*CritValCorrJohansenMultEJ[k,dets])
                UnionRejectionLocalAsyPowerEJAsymmMult[k,dets,cc,rr] <- mean(UnionRejectionStatEJAsymmMult)
                UnionRejectionStatBJAsymmMult <- (BoswijkStat*(BoswijkStat > CritvalBoswijk*CritValCorrBoswijkMultBJ[k,dets]) > CritvalBoswijk*CritValCorrBoswijkMultBJ[k,dets]) + (JohansenStat*(BoswijkStat < CritvalBoswijk*CritValCorrBoswijkMultBJ[k,dets]) > CritvalJohansen*CritValCorrJohansenMultBJ[k,dets])
                UnionRejectionLocalAsyPowerBJAsymmMult[k,dets,cc,rr] <- mean(UnionRejectionStatBJAsymmMult)
                UnionRejectionStatBERCAsymmMult <- (BoswijkStat*(BoswijkStat > CritvalBoswijk*CritValCorrBoswijkMultBERC[k,dets]) > CritvalBoswijk*CritValCorrBoswijkMultBERC[k,dets]) + (Re(ErrCorrStat)*(BoswijkStat < CritvalBoswijk*CritValCorrBoswijkMultBERC[k,dets]) < CritvalErrCorr*CritValCorrErrCorrMultBERC[k,dets])
                UnionRejectionLocalAsyPowerBERCAsymmMult[k,dets,cc,rr] <- mean(UnionRejectionStatBERCAsymmMult)

## -------------------------------- Fisher Type Tests --------------------------------"
                # Define Statistics" # ErrCorrPValue
                # 2 tests
                FisherStatEJ <- -2*(log(EngleGrangerPValue)+log(JohansenPValue))
                FisherStatBJ <- -2*(log(BoswijkPValue)+log(JohansenPValue))
                FisherStatBE <- -2*(log(BoswijkPValue)+log(EngleGrangerPValue))
                FisherStatBECR <- -2*(log(ErrCorrPValue)+log(BoswijkPValue))
                
                # 3 tests
                FisherStatBJE <- -2*(log(EngleGrangerPValue)+log(JohansenPValue)+log(BoswijkPValue))
                FisherStatBECRJ <- -2*(log(ErrCorrPValue)+log(BoswijkPValue)+log(JohansenPValue))
                FisherStatBECRE <- -2*(log(ErrCorrPValue)+log(BoswijkPValue)+log(EngleGrangerPValue)) # new Jens
                FisherStatECREJ <- -2*(log(ErrCorrPValue)+log(JohansenPValue)+log(EngleGrangerPValue)) # new Jens
                #4 tests
                FisherStatBECRJE <- -2*(log(ErrCorrPValue)+log(BoswijkPValue)+log(EngleGrangerPValue)+log(JohansenPValue))
                
                
                FisherStatEJw=-4*((1/(1+exp(5*(R2run-.25))))*log(EngleGrangerPValue)+(1-1/(1+exp(5*(R2run-.25))))*log(JohansenPValue))
                
                if (cc == 1 && rr == 1){# Write Critical Values"
                    # 2 tests 
                    NullDistrFisherEJ <- sort(FisherStatEJ)
                    CritvalFisherEJ <- NullDistrFisherEJ[rep*.95]
                    CritvalFisherEJ1 <- NullDistrFisherEJ[rep*.99]
                    CritvalFisherEJ2 <- NullDistrFisherEJ[rep*.95]
                    CritvalFisherEJ3 <- NullDistrFisherEJ[rep*.90]

                    NullDistrFisherBJ <- sort(FisherStatBJ)
                    CritvalFisherBJ <- NullDistrFisherBJ[rep*.95]
                    CritvalFisherBJ1 <- NullDistrFisherBJ[rep*.99]
                    CritvalFisherBJ2 <- NullDistrFisherBJ[rep*.95]
                    CritvalFisherBJ3 <- NullDistrFisherBJ[rep*.90]

                    NullDistrFisherBE <- sort(FisherStatBE)
                    CritvalFisherBE <- NullDistrFisherBE[rep*.95]
                    CritvalFisherBE1 <- NullDistrFisherBE[rep*.99]
                    CritvalFisherBE2 <- NullDistrFisherBE[rep*.95]
                    CritvalFisherBE3 <- NullDistrFisherBE[rep*.90]
                    
                    #3 Tests
                    NullDistrFisherBJE <- sort(FisherStatBJE)
                    CritvalFisherBJE <- NullDistrFisherBJE[rep*.95]
                    CritvalFisherBJE1 <- NullDistrFisherBJE[rep*.99]
                    CritvalFisherBJE2 <- NullDistrFisherBJE[rep*.95]
                    CritvalFisherBJE3 <- NullDistrFisherBJE[rep*.90]

                    NullDistrFisherBECR <- sort(FisherStatBECR)
                    CritvalFisherBECR <- NullDistrFisherBECR[rep*.95]
                    CritvalFisherBECR1 <- NullDistrFisherBECR[rep*.99]
                    CritvalFisherBECR2 <- NullDistrFisherBECR[rep*.95]
                    CritvalFisherBECR3 <- NullDistrFisherBECR[rep*.90]

                    NullDistrFisherBECRJE <- sort(FisherStatBECRJE)
                    CritvalFisherBECRJE <- NullDistrFisherBECRJE[rep*.95]
                    CritvalFisherBECRJE1 <- NullDistrFisherBECRJE[rep*.99]
                    CritvalFisherBECRJE2 <- NullDistrFisherBECRJE[rep*.95]
                    CritvalFisherBECRJE3 <- NullDistrFisherBECRJE[rep*.90]

                    NullDistrFisherBECRJ <- sort(FisherStatBECRJ)
                    CritvalFisherBECRJ <- NullDistrFisherBECRJ[rep*.95]
                    CritvalFisherBECRJ1 <- NullDistrFisherBECRJ[rep*.99]
                    CritvalFisherBECRJ2 <- NullDistrFisherBECRJ[rep*.95]
                    CritvalFisherBECRJ3 <- NullDistrFisherBECRJ[rep*.90]

                    NullDistrFisherEJw <- sort(FisherStatEJw)
                    CritvalFisherEJw <- NullDistrFisherEJw[rep*.95]
                    
                    ### Jens changes 
                    
                    #3 Tests
                    NullDistrFisherBECRE <- sort(FisherStatBECRE)
                    CritvalFisherBECRE  <- NullDistrFisherBECRE[rep*.95]
                    CritvalFisherBECRE1 <- NullDistrFisherBECRE[rep*.99]
                    CritvalFisherBECRE2 <- NullDistrFisherBECRE[rep*.95]
                    CritvalFisherBECRE3 <- NullDistrFisherBECRE[rep*.90]
                    
                    NullDistrFisherECREJ <- sort(FisherStatECREJ) 
                    CritvalFisherECREJ  <- NullDistrFisherECREJ[rep*.95]
                    CritvalFisherECREJ1 <- NullDistrFisherECREJ[rep*.99]
                    CritvalFisherECREJ2 <- NullDistrFisherECREJ[rep*.95]
                    CritvalFisherECREJ3 <- NullDistrFisherECREJ[rep*.90]
                    
                    
               }

                # Evaluate Local Power"
                FisherLocalAsyPowerEJ[k,dets,cc,rr] <- mean((FisherStatEJ > CritvalFisherEJ))
                FisherLocalAsyPowerBJ[k,dets,cc,rr] <- mean((FisherStatBJ > CritvalFisherBJ))
                FisherLocalAsyPowerBE[k,dets,cc,rr] <- mean((FisherStatBE > CritvalFisherBE))
                FisherLocalAsyPowerBJE[k,dets,cc,rr] <- mean((FisherStatBJE > CritvalFisherBJE))
                FisherLocalAsyPowerBECR[k,dets,cc,rr] <- mean((FisherStatBECR > CritvalFisherBECR))
                FisherLocalAsyPowerBECRJE[k,dets,cc,rr] <- mean((FisherStatBECRJE > CritvalFisherBECRJE))
                FisherLocalAsyPowerBECRJ[k,dets,cc,rr] <- mean((FisherStatBECRJ > CritvalFisherBECRJ))
                FisherLocalAsyPowerEJw[k,dets,cc,rr] <- mean((FisherStatEJw > CritvalFisherEJw))

                ## -------------------------------- Inverse Normal -------------------------------- */"
#                InvNormStatBJE <- 1/sqrt(3)*(norminv(EngleGrangerPValue)+norminv(JohansenPValue)+norminv(BoswijkPValue))"
#                if (cc == 1 && rr == 1){
#                    NullDistrInvNormBJE <- sort(InvNormStatBJE)
#                    CritvalInvNormBJE <- NullDistrInvNormBJE[rep*.05]
#                    CritvalInvNormBJE1 <- NullDistrInvNormBJE[rep*.01]
#                    CritvalInvNormBJE2 <- NullDistrInvNormBJE[rep*.05]
#                    CritvalInvNormBJE3 <- NullDistrInvNormBJE[rep*.10]
#                }
#                InvNormLocalAsyPowerBJE[k,dets,cc,rr]=mean((InvNormStatBJE <= CritvalInvNormBJE))  #condition not right
                ## -------------------------------- Min Pval -------------------------------- *"
                # Define Test Statistics"
                MinStatBECR <- min(ErrCorrPValue,BoswijkPValue)
                MinStatEJ <- min(EngleGrangerPValue,JohansenPValue)

                if (cc == 1 && rr == 1){# Calculate Critical Values"
                    NullDistrMinBECR <- sort(MinStatBECR)
                    CritvalMinBECR <- NullDistrMinBECR[rep*.05]
                    CritvalMinBECR1 <- NullDistrMinBECR[rep*.01]
                    CritvalMinBECR2 <- NullDistrMinBECR[rep*.05]
                    CritvalMinBECR3 <- NullDistrMinBECR[rep*.10]

                    NullDistrMinEJ <- sort(MinStatEJ)
                    CritvalMinEJ <- NullDistrMinEJ[rep*.05]
                    CritvalMinEJ1 <- NullDistrMinEJ[rep*.01]
                    CritvalMinEJ2 <- NullDistrMinEJ[rep*.05]
                    CritvalMinEJ3 <- NullDistrMinEJ[rep*.10]
                }

                # Local Power"
                MinLocalAsyPowerBERC[k,dets,cc,rr] <- mean((MinStatBECR < CritvalMinBECR))
                MinLocalAsyPowerEJ[k,dets,cc,rr] <- mean((MinStatEJ < CritvalMinEJ))

            }# R2"
        }# c"

        # -------------------------------- results -------------------------------- */"
#        disp (['k && case ' num2str(k) num2str(dets)])
#        disp (['c.v. Boswijk ' num2str(CritvalBoswijk)])
#        disp (['c.v. Johansen ' num2str(CritvalJohansen)])
#        disp (['c.v. EngleGranger ' num2str(CritvalEngleGranger)])
#        disp (['c.v. ErrCorr ' num2str(CritValErrCorr1)])
#        disp ' '"
#        disp (['c.v. FisherEJ ' num2str(CritvalFisherEJ)])
#        disp (['c.v. FisherBJ ' num2str(CritvalFisherBJ)])
#        disp (['c.v. FisherBE ' num2str(CritvalFisherBE)])
#        disp (['c.v. FisherBJE ' num2str(CritvalFisherBJE)])
#        disp (['c.v. FisherBECR ' num2str(CritvalFisherBECR)])
#        disp (['c.v. FisherEJ wght ' num2str(CritvalFisherEJw)])
#        disp (['c.v. InvNormBJE ' num2str(CritvalInvNormBJE)])
#        disp ' '
#        disp (['the chi(2*2) 95% quantile is' num2str(chi2inv(.95,4))])
#        disp (['the chi(2*3) 95% quantile is' num2str(chi2inv(.95,6))])
#        disp ' '
#        CV[,k,dets,1] <- CritvalBoswijk1; CritValErrCorr1; CritvalJohansen1; CritvalEngleGranger1; CritvalFisherBECR1
#            CritvalFisherBJ1; CritvalFisherBE1; CritvalFisherEJ1; CritvalFisherBJE1; CritvalFisherBECRJ1; CritvalFisherBECRJE1; CritvalInvNormBJE1; CritvalMinBECR1; CritvalMinEJ1]
#        CV(:,k,dets,2) <- [CritvalBoswijk2; CritValErrCorr2; CritvalJohansen2; CritvalEngleGranger2; CritvalFisherBECR2; ...
#            CritvalFisherBJ2; CritvalFisherBE2; CritvalFisherEJ2; CritvalFisherBJE2; CritvalFisherBECRJ2; CritvalFisherBECRJE2; CritvalInvNormBJE2; CritvalMinBECR2; CritvalMinEJ2]
#        CV(:,k,dets,3) <- [CritvalBoswijk3; CritValErrCorr3; CritvalJohansen3; CritvalEngleGranger3; CritvalFisherBECR3; ...
#            CritvalFisherBJ3; CritvalFisherBE3; CritvalFisherEJ3; CritvalFisherBJE3; CritvalFisherBECRJ3; CritvalFisherBECRJE3; CritvalInvNormBJE3; CritvalMinBECR3; CritvalMinEJ3]
#        NullDistr[ , ,k , dets] <- cbind(NullDistrEngleGranger, NullDistrJohansen, NullDistrErrCorr, NullDistrBoswijk)

    }# dets */"
}# k */"
#save workspace
#if (savecv==1){
#   save critical_values CV
#    save NullDistributions NullDistr"
#}


##### Am Ende muss du alles in einer RDATA datei laden 


