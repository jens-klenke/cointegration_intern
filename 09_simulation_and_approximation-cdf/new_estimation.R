#---- Preliminary ---- 
# packages
source(here::here('01_code/packages/packages.R'))

# load functions 
source(here::here('09_simulation_and_approximation-cdf/estimation_functions.R'))

#-- tibble with models and data ----
## Load Simulation Data 
if(Sys.info()['nodename'] == "DELL-ARBEIT") { # Jens 
    Data <- readRDS('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_100k.rds')
    # load('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_1_m.RData')
} else if(Sys.info()['user'] == "Janine") { # Janine
    load("/Users/Janine/Desktop/BayerHanck/Data_1_m.RData")
} else if(Sys.info()['nodename'] == "OEK-TS01") { # Server
    load('D:\\Klenke\\Data_1_m.RData')
}

# just for the programming, delete before running final results
 Data %<>%
    dplyr::sample_n(1000000)

# Split Dataset in Cases
data_case_1 <- Data %>%
    dplyr::filter(case == 1)
data_case_2 <- Data %>%
    dplyr::filter(case == 2)
data_case_3 <- Data %>%
    dplyr::filter(case == 3)

rm('Data')
#---- Boxcox Transformation ----

# case_1 
data_case_1 %<>% bc_log_E_J_fun()

# case_2 
data_case_2 %<>% bc_log_E_J_fun()

# case_3 
data_case_3 %<>% bc_log_E_J_fun()

#lambda_bc_EJ <- tibble(
#    case = rep(1:3, each = 2), 
#    side = rep(c('stat', 'p'), times = 3),
#    value = c(lambda_stat_case_1, lambda_p, lambda_stat_case_2, lambda_p,  lambda_stat_case_3, lambda_p)   
#)

#-- calls and range of power ----
# E_J
calls_E_J <- c('p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, power)',
#           'p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, power) + I(1/k)',
#           'p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, power) + k + I(1/k)',
           'p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, power) + I(1/k) + poly(stat_Fisher_E_J, power) * I(1/k)',
#           'p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, power) + log(k) + poly(stat_Fisher_E_J, power) * log(k)',
#           'p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), power) + k',
#           'p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), power) * k',
#           'p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), power) * k + (1/k)',
#           'p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), power) + log(k)',
           'p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), power) * log(k)'
           )

# all
calls_all <- c('p_value_Fisher_all ~ poly(stat_Fisher_all, power)')

# power
expo <- 3:5

#--------------- E_J -----
#-- E_J case_1 ----
table_E_J_case_1 <- table_E_J_fun(data_case_1)
# bei implementierung noch die plot values löschen!

#-- E_J case_2 ----
table_E_J_case_2 <- table_E_J_fun(data_case_2)

#-- E_J case_3 ----
table_E_J_case_3 <- table_E_J_fun(data_case_3)


#--------------- ALL -----
#-- all case_1 ----
table_all_case_1 <- table_all_fun(data_case_1)

#-- all case_2 ----
table_all_case_2 <- table_all_fun(data_case_2)

#-- all case_3 ----
table_all_case_3 <- table_all_fun(data_case_3)


