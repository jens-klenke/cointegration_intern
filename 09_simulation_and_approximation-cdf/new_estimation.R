#---- Preliminary ---- 
# packages
source(here::here('01_code/packages/packages.R'))

# load functions 
source(here::here('09_simulation_and_approximation-cdf/estimation_functions.R'))

# parallel
plan(multisession, workers = 2)
options(future.globals.maxSize = 2.147e+9)

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
    dplyr::sample_n(100000)

# Split Dataset in Cases
data_case_1 <- Data %>%
    dplyr::filter(case == 1)
data_case_2 <- Data %>%
    dplyr::filter(case == 2)
data_case_3 <- Data %>%
    dplyr::filter(case == 3)

rm('Data')
#---- Boxcox Transformation ----

# Hier bitte kein Piping!

# case_1 
data_case_1 <- bc_log_fun(data_case_1)

# case_2 
data_case_2 <- bc_log_fun(data_case_2)

# case_3 
data_case_3 <- bc_log_fun(data_case_3)

lambda_bc_EJ <- tibble(
    case = rep(1:3, each = 2), 
    side = rep(c('stat', 'p'), times = 3),
    value = c(lambda_stat_E_J_1, lambda_p, lambda_stat_E_J_2, lambda_p, 
              lambda_stat_E_J_3, lambda_p)   
)

lambda_bc_all <- tibble(
    case = rep(1:3, each = 2), 
    side = rep(c('stat', 'p'), times = 3),
    value = c(lambda_stat_all_1, lambda_p, lambda_stat_all_2, lambda_p, 
              lambda_stat_all_3, lambda_p)   
)

#-- calls and range of power ----
# E_J 
calls_E_J <- c('p_value_Fisher ~ poly(stat_Fisher_E_J, power)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J, power) + k',
               'p_value_Fisher ~ poly(stat_Fisher_E_J, power) * k' ,
               'p_value_Fisher ~ poly(stat_Fisher_E_J, power) * k * I(1/k)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J, power) * log(k)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J, power) * log(k) * I(1/k)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J, power) * log(k) * I(1/k) * k',
               'p_value_Fisher ~ poly(stat_Fisher_E_J, power) * log(k) * I(1/k) * k * sqrt(k)',
               # bc
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) + k',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) * k',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) * k * I(1/k)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) * log(k)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) * log(k) * I(1/k)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) * log(k) * I(1/k) * k',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) * log(k) * I(1/k) * k * sqrt(k)',
               # lg
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J, power)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J, power) + k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J, power) * k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J, power) * k * I(1/k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J, power) * log(k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J, power) * log(k) * I(1/k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J, power) * log(k) * I(1/k) * k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J, power) * log(k) * I(1/k) * k * sqrt(k)',
               # lg + bc
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) + k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) * k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) * k * I(1/k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) * log(k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) * log(k) * I(1/k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) * log(k) * I(1/k) * k', 
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) * log(k) * I(1/k) * k * sqrt(k)', 
               # bc + bc
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) + k',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) * k',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) * k * I(1/k)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) * log(k)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) * log(k) * I(1/k)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) * log(k) * I(1/k) * k', 
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) * log(k) * I(1/k) * k * sqrt(k)'
           )

# all
calls_all <- c('p_value_Fisher ~ poly(stat_Fisher_all, power)',
               'p_value_Fisher ~ poly(stat_Fisher_all, power) + k',
               'p_value_Fisher ~ poly(stat_Fisher_all, power) * k',
               'p_value_Fisher ~ poly(stat_Fisher_all, power) * k * I(1/k)',
               'p_value_Fisher ~ poly(stat_Fisher_all, power) * log(k)',
               'p_value_Fisher ~ poly(stat_Fisher_all, power) * log(k) * I(1/k)',
               'p_value_Fisher ~ poly(stat_Fisher_all, power) * log(k) * I(1/k) * k',
               'p_value_Fisher ~ poly(stat_Fisher_all, power) * log(k) * I(1/k) * k * sqrt(k)',
               # bc
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power)',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) + k',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) * k',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) * k * I(1/k)',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) * log(k)',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) * log(k) * I(1/k)',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) * log(k) * I(1/k) * k',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) * log(k) * I(1/k) * k * sqrt(k)',
               # lg
               'p_value_Fisher_lg ~ poly(stat_Fisher_all, power)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all, power) + k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all, power) * k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all, power) * k * I(1/k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all, power) * log(k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all, power) * log(k) * I(1/k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all, power) * log(k) * I(1/k) * k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all, power) * log(k) * I(1/k) * k * sqrt(k)',
               # lg + bc
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) + k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) * k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) * k * I(1/k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) * log(k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) * log(k) * I(1/k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) * log(k) * I(1/k) * k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) * log(k) * I(1/k) * k * sqrt(k)', 
               # bc + bc
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) + k',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) * k',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) * k * I(1/k)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) * log(k)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) * log(k) * I(1/k)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) * log(k) * I(1/k) * k', 
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) * log(k) * I(1/k) * k * sqrt(k)'
)

# power
expo <- 3:10

#--------------- E_J -----
#-- E_J case_1 ----
tictoc::tic()
table_E_J_case_1 <- table_E_J_fun(data_case_1)
tictoc::toc()

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


