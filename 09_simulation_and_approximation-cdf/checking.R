#---- Preliminary ---- 
# packages
source(here::here('01_code/packages/packages.R'))

#-- data ----
## Load Simulation Data 

# p-values vs cdf   
cowplot::plot_grid(
Data %>%
    dplyr::filter(case == 1,
                  k == 1) %>%
    ggplot2::ggplot(aes(x = stat_Fisher_all, y = p_value_Fisher_all)) +
    ggplot2::geom_line(),

Data %>%
    dplyr::filter(case == 1,
                  k == 1) %>%
    dplyr::mutate(cdf = seq(1/100000,1, by = 1/100000)) %>%
    ggplot2::ggplot(aes(x = stat_Fisher_all, y = cdf)) +
    ggplot2::geom_line())

# wahrscheinlich klammer fehler in der SIM! 

Fisher_all_critical_val <- Data %>%
    dplyr::group_by(case, k) %>%
    dplyr::filter(p_value_Fisher_all == 0.9500) %>%
    dplyr::select(stat_Fisher_all, p_value_Fisher_all, k, case)


Fisher_E_J_critical_val <- Data %>%
    dplyr::group_by(case, k) %>%
    dplyr::filter(p_value_Fisher_E_J == 0.9500) %>%
    dplyr::select(stat_Fisher_E_J, p_value_Fisher_E_J, k, case)

rm('Data')
# load values from the package
load("~/GitHub/bayerhanck/R/sysdata.rda")

## E_J 
bh <- as_tibble(crit_val_1_0.05) 

bh %<>% 
    dplyr::rename('1' = V1,
                  '2' = V2,
                  '3' = V3) %>%
    dplyr::mutate(k = 1:11) %>%
    tidyr::pivot_longer(col = c('1', '2', '3'), names_to = 'case') %>%
    dplyr::mutate(case = as.integer(case))

compare_E_J <- Fisher_E_J_critical_val %>%
    left_join( bh, by = c('case', 'k')) %>%
    dplyr::mutate(diff = stat_Fisher_E_J - value)

## all
bh_2 <- as_tibble(crit_val_2_0.05) 

bh_2 %<>% 
    dplyr::rename('1' = V1,
                  '2' = V2,
                  '3' = V3) %>%
    dplyr::mutate(k = 1:11) %>%
    tidyr::pivot_longer(col = c('1', '2', '3'), names_to = 'case') %>%
    dplyr::mutate(case = as.integer(case))

compare_all <- Fisher_all_critical_val %>%
    left_join( bh_2, by = c('case', 'k')) %>%
    dplyr::mutate(diff = stat_Fisher_all - value)


#### data wrangling 


#---------------------------------------------------- Janine  -----------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#-- tibble with models and data ----
## Load Simulation Data 
if(Sys.info()['nodename'] == "DELL-ARBEIT") { # Jens 
    Data <- readRDS('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_100k.rds')
    # load('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_1_m.RData')
} else if(Sys.info()['nodename'] == "MacBook-Pro.local") { # Janine
    load("/Users/Janine/Desktop/BayerHanck/Data_1_m.RData")
} else if(Sys.info()['nodename'] == "OEK-TS01") { # Server
    load('D:\\Klenke\\Data_1_m.RData')
}


Data %<>%
    dplyr::rename('cdf_Fisher_all' = p_value_Fisher_all,
                  'cdf_Fisher_E_J' = p_value_Fisher_E_J) %>%
    dplyr::mutate(p_value_Fisher_all = 1 - cdf_Fisher_all,
                  p_value_Fisher_E_J = 1 - cdf_Fisher_E_J) %>%
    dplyr::select(-c(cdf_Fisher_all, cdf_Fisher_E_J))

save(Data, file = "/Users/Janine/Desktop/BayerHanck/Data_1_m.RData")


