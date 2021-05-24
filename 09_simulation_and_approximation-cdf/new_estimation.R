#---- Preliminary ---- 
# packages
source(here::here('01_code/packages/packages.R'))

# load functions 
source(here::here('09_simulation_and_approximation-cdf/estimation_functions.R'))

# parallel
#doParallel::registerDoParallel(cores = 6)
plan(multisession, workers = 3)
# set off warning for random numbers
options(future.rng.onMisuse = "ignore")

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

# k as factor/dummy
Data %<>% dplyr::mutate(
    k_dummy = as.factor(k)
    )

# just for programming, delete before running final results
# Data %<>%
#      dplyr::sample_n(3300000)

# Split Dataset in Cases
data_case_1 <- Data %>%
    dplyr::filter(case == 1)
data_case_2 <- Data %>%
    dplyr::filter(case == 2)
data_case_3 <- Data %>%
    dplyr::filter(case == 3)

rm('Data')

#---- Boxcox Transformation ----
data_case_1 <- bc_log_fun(data_case_1)
data_case_2 <- bc_log_fun(data_case_2)
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

lambda_values <- bind_rows(lambda_bc_all, 
                           lambda_bc_EJ) %>% 
    dplyr::mutate(test = c(rep("all", 6), 
                           rep("e_j", 6)))

# save(lambda_values, file = here::here('09_simulation_and_approximation-cdf/lambda_package.RData'))


#---- calls and range of power ----
# E_J 
calls_E_J <- c('p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) + k',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) * k',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) + log(k)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) * log(k)',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) + k_dummy',
               'p_value_Fisher ~ poly(stat_Fisher_E_J_bc, power) * k_dummy',
               # lg
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) + k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) * k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) + log(k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) * log(k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) + k_dummy',
               'p_value_Fisher_lg ~ poly(stat_Fisher_E_J_bc, power) * k_dummy',
               # bc
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) + k',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) * k',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) + log(k)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) * log(k)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) + k_dummy',
               'p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, power) * k_dummy'
               )

# all
calls_all <- c('p_value_Fisher ~ poly(stat_Fisher_all_bc, power)',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) + k',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) * k',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) + log(k)',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) * log(k)',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) + k_dummy',
               'p_value_Fisher ~ poly(stat_Fisher_all_bc, power) * k_dummy',
               # lg
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) + k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) * k',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) + log(k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) * log(k)',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) + k_dummy',
               'p_value_Fisher_lg ~ poly(stat_Fisher_all_bc, power) * k_dummy',
               # bc
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) + k',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) * k',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) + log(k)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) * log(k)',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) + k_dummy',
               'p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, power) * k_dummy'
               )

# power
expo <- 12:13


#---- E_J ----
table_E_J_case_1 <- table_fun(data_case_1, "E_J")
table_E_J_case_2 <- table_fun(data_case_2, "E_J")
table_E_J_case_3 <- table_fun(data_case_3, "E_J")

#---- all ----
table_all_case_1 <- table_fun(data_case_1, "all")
table_all_case_2 <- table_fun(data_case_2, "all")
table_all_case_3 <- table_fun(data_case_3, "all")

save(table_E_J_case_1, table_E_J_case_2, table_E_J_case_3,
     table_all_case_1, table_all_case_2, table_all_case_3,
     file = here::here('09_simulation_and_approximation-cdf/server_results.RData'))
