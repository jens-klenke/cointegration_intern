#---- Preliminary ---- 
# packages 
source(here::here('01_code/packages/packages.R'))

# functions 
source(here::here('09_simulation_and_approximation-cdf/func_p_approximation.R'))

# load models 
load(here::here('09_simulation_and_approximation-cdf/models_package.RData'))

# lambda for boc-cox
load(here::here('09_simulation_and_approximation-cdf/lambda_package.RData'))

# vom test selbst geliefert 
k = 2 
bh.test <- a$bh.test
trendtype  <- 1
test.type <- a$test.type


get_p_value <- function(bh.test, trendtype, test.type){

#lambda_p <<- get_lambda(trendtype, 'p', 'all')


new_data <- tibble(stat_Fisher_all = bh.test, 
                   stat_Fisher_E_J = bh.test, 
                   k =  k )

new_data %<>%
    dplyr::mutate(stat_Fisher_E_J_bc = ((stat_Fisher_E_J^get_lambda(trendtype, 'stat', 'e_j'))-1)/get_lambda(trendtype, 'stat', 'e_j'),
                  stat_Fisher_all_bc = ((stat_Fisher_all^get_lambda(trendtype, 'stat', 'all'))-1)/get_lambda(trendtype, 'stat', 'all'))

    

p.value_raw <- suppressWarnings(predict(get_model(trendtype, test.type), new_data))

p.value_trans <- if(get_p_trans(trendtype, test.type) == 'log'){exp(p.value_raw)} else{invBoxCox(p.value_raw)}

p.value <- ifelse(p.value_trans >= 1, 9.9999e-1, ifelse(p.value_trans <= 0, 1e-12, p.value_trans))

return(p.value)

}

get_p_value(20.1, 1, 'all')
get_p_value(50, 1, 'all')
get_p_value(55, 1, 'all')
get_p_value(100, 1, 'all')
get_p_value(200, 1, 'all')



## Janine ausführen und speichern 

# Load Simulation Data 
if(Sys.info()['nodename'] == "DELL-ARBEIT") { # Jens 
    Data <- readRDS('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_100k.rds')
    # load('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_1_m.RData')
} else if(Sys.info()['user'] == "Janine") { # Janine
    load("/Users/Janine/Desktop/BayerHanck/Data_1_m.RData")
} else if(Sys.info()['nodename'] == "OEK-TS01") { # Server
    load('D:\\Klenke\\Data_1_m.RData')
}

crit_val_all <- Data %>%
    dplyr::group_by(case, k) %>%
    dplyr::select(stat_Fisher_all) %>%
    dplyr::slice_max(stat_Fisher_all , n = 100) %>%
    dplyr::filter(stat_Fisher_all == min(stat_Fisher_all))

crit_val_e_j <- Data %>%
    dplyr::group_by(case, k) %>%
    dplyr::select(stat_Fisher_E_J) %>%
    dplyr::slice_max(stat_Fisher_E_J , n = 100) %>%
    dplyr::filter(stat_Fisher_E_J == min(stat_Fisher_E_J))

