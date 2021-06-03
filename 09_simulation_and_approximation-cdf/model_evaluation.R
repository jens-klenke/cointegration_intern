#---- Preliminary ---- 
# packages  
source(here::here('01_code/packages/packages.R'))

# functions 
source(here::here('09_simulation_and_approximation-cdf/functions_model_eval.R'))

# load metrics
load(here::here('09_simulation_and_approximation-cdf/server_results.RData'))

# final models and lambda
load(here::here('09_simulation_and_approximation-cdf/lambda_package.RData'))

# Load Data
if(Sys.info()['nodename'] != "DELL-ARBEIT") {
    load(here::here('09_simulation_and_approximation-cdf/data_cases.RData'))
} else if(Sys.info()['nodename'] == "DELL-ARBEIT") { # Jens 
    load('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\data_cases.RData')
}


lambda_p <- get_lambda(lambda_values, 1, 'p', 'all')


#---- Analysis tables ----
sum_table_all_case_1 <- best_5_table(table_all_case_1)

sum_table_all_case_2 <- best_5_table(table_all_case_2)

sum_table_all_case_3 <- best_5_table(table_all_case_3)

sum_table_E_J_case_1 <- best_5_table(table_E_J_case_1)

sum_table_E_J_case_2 <- best_5_table(table_E_J_case_2)

sum_table_E_J_case_3 <- best_5_table(table_E_J_case_3)


# ---- Plots ----
# Cut Data 
data_case_1_small <- data_case_1 %>%
    dplyr::filter(p_value_Fisher %in% c(min(p_value_Fisher), seq(0, 1, 0.0001)))
data_case_2_small <- data_case_2  %>% 
    dplyr::filter(p_value_Fisher %in% c(min(p_value_Fisher), seq(0, 1, 0.0001)))
data_case_3_small <- data_case_3 %>% 
    dplyr::filter(p_value_Fisher %in% c(min(p_value_Fisher), seq(0, 1, 0.0001)))


# Add Predictions
# case = 1
data_case_1_all <- add_pred(data_case_1_small, 'all')
data_case_1_E_J <- add_pred(data_case_1_small, 'E_J')

# case = 2
data_case_2_all <- add_pred(data_case_2_small, 'all')
data_case_2_E_J <- add_pred(data_case_2_small, 'E_J')

# case = 3
data_case_3_all <- add_pred(data_case_3_small, 'all')
data_case_3_E_J <- add_pred(data_case_3_small, 'E_J')

# Plot approximated against simulated p-values
# case = 1
data_case_1_all %>% own_plot()
data_case_1_E_J %>% own_plot()

# case = 2
data_case_2_all %>% own_plot()
data_case_2_E_J %>% own_plot()

# case = 3
data_case_3_all %>% own_plot()
data_case_3_E_J %>% own_plot()


# Plot approximated against simulated p-values <= 0.2
# case = 1
data_case_1_all %>% own_plot_0.2()
data_case_1_E_J %>% own_plot_0.2()

# case = 2
data_case_2_all %>% own_plot_0.2()
data_case_2_E_J %>% own_plot_0.2()

# case = 3
data_case_3_all %>% own_plot_0.2()
data_case_3_E_J %>% own_plot_0.2()

