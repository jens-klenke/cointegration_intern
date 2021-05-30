#---- Preliminary ---- 
# packages 
source(here::here('01_code/packages/packages.R'))

# functions 
source(here::here('09_simulation_and_approximation-cdf/functions_model_eval.R'))

# load metrics
load(here::here('09_simulation_and_approximation-cdf/server_results.RData'))

# final models and lambda
load(here::here('09_simulation_and_approximation-cdf/lambda_package.RData'))


#models <- tibble(
#    test = rep(c('all', 'e_j'), times = c(3, 3)),
#    case = rep(1:3, times = 2),
#    models = list(clean_mod_all_1, clean_mod_all_2, clean_mod_all_3,
#                  clean_mod_E_J_1, clean_mod_E_J_2, clean_mod_E_J_3), 
#    response = c('bc', 'bc', 'bc', 'log', 'log', 'log')
#)

#save(models, file = here::here('09_simulation_and_approximation-cdf/models_package.RData'))

# Load Simulation Data 
if(Sys.info()['nodename'] == "DELL-ARBEIT") { # Jens 
    Data <- readRDS('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_100k.rds')
    # load('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_1_m.RData')
} else if(Sys.info()['user'] == "Janine") { # Janine
    load("/Users/Janine/Desktop/BayerHanck/Data_1_m.RData")
} else if(Sys.info()['nodename'] == "OEK-TS01") { # Server
    load('D:\\Klenke\\Data_1_m.RData')
}

lambda_p <- get_lambda(lambda_values, 1, 'p', 'all')

#-- Analysis  tables ----

sum_table_all_case_1 <- best_5_table(table_all_case_1)

sum_table_all_case_2 <- best_5_table(table_all_case_2)

sum_table_all_case_3 <- best_5_table(table_all_case_3)

sum_table_E_J_case_1 <- best_5_table(table_E_J_case_1)

sum_table_E_J_case_2 <- best_5_table(table_E_J_case_2)

sum_table_E_J_case_3 <- best_5_table(table_E_J_case_3)

# Cut Data 
Data %<>%
    dplyr::select(-p_value_Fisher_all) %>%
    dplyr::mutate(p_value_Fisher_bc = ((p_value_Fisher_E_J^lambda_p)-1)/lambda_p, 
                  p_value_Fisher_lg = log(p_value_Fisher_E_J)) %>%
    dplyr::rename(p_value_Fisher = p_value_Fisher_E_J) %>%
    dplyr::filter(p_value_Fisher %in% c(min(p_value_Fisher), seq(0, 1, 0.0001))) 

### 

# plot data
# case = 1
data_case_1_all <- plot_data(Data, 1, 'all')
data_case_1_e_j <- plot_data(Data, 1, 'e_j')

# case = 2
data_case_2_all <- plot_data(Data, 2, 'all')
data_case_2_e_j <- plot_data(Data, 2, 'e_j')

# case = 3
data_case_3_all <- plot_data(Data, 3, 'all')
data_case_3_e_j <- plot_data(Data, 3, 'e_j')




#-- Analysis  plots ----

### case 1

# all
own_plot(data_case_1_all)

# anderes modell
data_case_1_all %>%
    dplyr::filter(p_value_Fisher <= 0.2) %>%
    own_plot(max_graph = 0.2)



# e_j
own_plot(data_case_1_e_j)

data_case_2_e_j %>%
    dplyr::filter(p_value_Fisher <= 0.2) %>%
    own_plot(max_graph = 0.2)

### case 2

# all
own_plot(data_case_2_all)

# eventuell anderes Modell 
data_case_2_all %>%
    dplyr::filter(p_value_Fisher <= 0.2) %>%
    own_plot(max_graph = 0.2)

# e_j
own_plot(data_case_2_e_j)

data_case_2_e_j %>%
    dplyr::filter(p_value_Fisher <= 0.2) %>%
    own_plot(max_graph = 0.2)

### case 3

# all
own_plot(data_case_3_all)

data_case_3_all %>%
    dplyr::filter(p_value_Fisher <= 0.2) %>%
    own_plot(max_graph = 0.2)

# e_j
own_plot(data_case_3_e_j)

data_case_3_e_j %>%
    dplyr::filter(p_value_Fisher <= 0.2) %>%
    own_plot(max_graph = 0.2)


