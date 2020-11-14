## packages
source(here::here('01_code/packages/packages.R'))

## sub functions 
source(here::here('09_simulation_and_approximation-cdf/sim_functions.R'))

## parallelation 
num_cores <- detectCores() # need to change for the server
use_cores <- num_cores-1  # need to change for the server

registerDoParallel(use_cores) # cores need to be hard coded for the server

## Parameters
R2 <- seq(0, 0.95, 0.05) # long term correlation
N <- 1000 # length for every trajectory
c_run <- 0 # Parameter pesavanto model
K <- 11 #  max of Variables
rep <- 1000000 # Number of Repetitions 100000
cases <- 3 # cases 
lambda <- (seq(1:N)/N) 

## Loop
tictoc::tic() 
#for (k in 1:K){# Number of Regressor Loop ### parallelisieren
Data <- foreach (k = 1:K, .combine = rbind) %dopar% { 
    data <- NULL   
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
        BoswijkPValu <- rep(NA, rep)
        JohansenPValue <- rep(NA, rep)
        EngleGrangerPValue <- rep(NA, rep)
        ErrCorrPValue <- rep(NA, rep)
        
        for (rr in 1:length(R2)){ # Loop over Pesaventos R2"

                # Set R2 
                R2run <- R2[rr]
                
                # Loop over repetitions"
                for (j in 1:rep){
                    
                    # -----------------------------Common Terms------------#
                    
                    c_terms <- sim_terms(N, k, R2run, c_run, dets)
                    
                    # -------------------------------- Boswijk------------------------ #
                    NullStatBoswijk[j] <- stat_Boswijk(c_run, N, c_terms$J12dc_sq, c_terms$J12DW2, c_terms$WdcDW2, c_terms$WdcWdci)
                    
                    # -------------------------------- Johansen -------------------------------- #
                    
                    NullStatJohansen[j] <- stat_Johansen(Wdc = c_terms$Wdc, J12dc = c_terms$J12dc, u = c_terms$u, 
                                                         N = N, k = k, c_run = c_run, WdcWdci = c_terms$WdcWdci)
                    
                    
                    # -------------------------------- Engle-Granger -------------------------------- #
                    
                    NullStatEngleGranger[j] <- stat_Eng_Gr(c_run = c_run, R2run = R2run, N = N, k = k, W1dW1di = c_terms$W1dW1di, 
                                                           W1d = c_terms$W1d, J12dc = c_terms$J12dc, Wdc = c_terms$Wdc, 
                                                           u = c_terms$u, u12 = c_terms$u12)
        
                    # -------------------------------- ECR (Banerjee) -------------------------------- #
                
                    NullStatErrCorr[j] <- stat_Banerjee(c_run = c_run, R2run = R2run, N = N, k = k, J12DW2 = c_terms$J12DW2, 
                                                        W1dJ12dc = c_terms$W1dJ12dc, W1dW1di = c_terms$W1dW1di, 
                                                        W1d = c_terms$W1d, J12dc_sq = c_terms$J12dc_sq, u = c_terms$u)
                    
                }#rep loop end
                
                ## Write Null Distributions and p-values for the underlying tests 
                ## p-values are ranked by the NullStat and not by NULLDist to not lose the association
                NullDistrBoswijk <- sort(NullStatBoswijk)
                BoswijkPValue <- 1 - rank(NullStatBoswijk)/rep+10^(-100)
                
                NullDistrJohansen <- sort(NullStatJohansen)
                JohansenPValue <- 1- rank(NullStatJohansen)/rep+10^(-100)
                
                NullDistrEngleGranger <- sort(NullStatEngleGranger)
                EngleGrangerPValue <- rank(NullStatEngleGranger)/rep+10^(-100)
                
                NullDistrErrCorr <- sort(NullStatErrCorr)
                ErrCorrPValue <- rank(NullStatErrCorr)/rep+10^(-100)
                
                ## -------------------------------- Fisher Type Tests --------------------------------"
                
                # Fisher test statistics 
                # 2 tests 
                FisherStat_E_J <- -2*(log(EngleGrangerPValue) + log(JohansenPValue))
                
                #4 tests
                FisherStat_all <- -2*( log(ErrCorrPValue) + log(BoswijkPValue) + 
                                       log(EngleGrangerPValue) + log(JohansenPValue))
                
                # Write Null Distributions and p-values 
                    
                # Engel - Johansen
                NullDistrFisher_E_J <- sort(FisherStat_E_J)
                Fisher_E_J_PValue <- rankindx(NullDistrFisher_E_J, 1)/rep+10^(-1000)

                #4 tests
                NullDistrFisher_all <- sort(FisherStat_all)
                Fisher_all_PValue <-rankindx(NullDistrFisher_all, 1)/rep+10^(-1000)
                
                # sorting of the p-values of the underlying tests 
                BoswijkPValue <- sort(BoswijkPValue, decreasing = TRUE)
                JohansenPValue <- sort(JohansenPValue, decreasing = TRUE)
                EngleGrangerPValue <- sort(EngleGrangerPValue)
                ErrCorrPValue <- sort(ErrCorrPValue) 
                
            }# R2
        
        D <- tibble::tibble(
            p_value_E_G = EngleGrangerPValue,
            stat_E_G = NullDistrEngleGranger,
            p_value_J = JohansenPValue,
            stat_J = NullDistrJohansen,
            p_value_Bo = BoswijkPValue,
            stat_Bo = NullDistrBoswijk,
            p_value_Ba = ErrCorrPValue,
            stat_Ba = NullDistrErrCorr,
            p_value_Fisher_E_J = Fisher_E_J_PValue,
            stat_Fisher_E_J = NullDistrFisher_E_J,
            p_value_Fisher_all = Fisher_all_PValue,
            stat_Fisher_all = NullDistrFisher_all,
            k = k,
            case = dets
            
        )
        data <- rbind(data, D)
         
    }# dets 
    data
}# k 
tictoc::toc()

base::saveRDS(Data, file = here::here('09_simulation_and_approximation-cdf/Data_1m.rds'))

