#---- Preliminary ---- 
# packages 
source(here::here('01_code/packages/packages.R'))

### functions
best_5_table <- function(data, col){
    formattable(data%>%
                               dplyr::slice_min(RMSE_cor_0.2, n = 5)%>%
                               dplyr::select(-c(all_of(col), expo, model)),
                           list(
                               RMSE = color_tile("green", "red"),
                               RMSE_cor = color_tile("green", "red"),
                               RMSE_0.2 = color_tile("green", "red"),
                               RMSE_cor_0.2 = color_tile("green", "red")
                           ))
}

best_model <- function(data){
    data %<>% 
        dplyr::filter(RMSE_cor_0.2 == min(RMSE_cor_0.2)) %>%
        dplyr::select(formula)
}

invBoxCox <- function(x){
    x <- if (lambda_p == 0) exp(as.complex(x)) else (lambda_p*as.complex(x) + 1)^(1/lambda_p)
    return(Re(x))
}

# getting lambda
get_lambda <- function(data, case_w, art){
    data %>%
        dplyr::filter(case == case_w, 
                      side == art) %>%
        dplyr::select(value) %>%
        dplyr::pull()
}

own_plot <- function(data, max_graph = 1){
    
    data %>%
        ggplot(aes(x = p_value_Fisher, y = PRED_cor_all)) +
        geom_line(color = '#004c93') +
        geom_segment(aes(x = 0, xend = max_graph, y = 0, yend = max_graph), 
                     linetype = 'dashed', size = 1, color = 'grey') +
        xlim(c(0, max_graph))+
        ylim(c(0, max_graph))+
        labs(x = '\n Simulated p-values', y = 'Approximated p-values \n')+
        theme_bw()+
        facet_wrap(~k)+
        theme(panel.spacing = unit(1, "lines"),
              strip.background = element_rect(colour = 'black',
                                              fill = '#004c93'),
              strip.text.x = element_text(size = 12, color = 'white' # , face = "bold.italic"
              ),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20))
}

# load metrics
load(here::here('09_simulation_and_approximation-cdf/server_results.RData'))

# final models and lambda
load(here::here('09_simulation_and_approximation-cdf/lambda_bc_ALL_server.RData'))
load(here::here('09_simulation_and_approximation-cdf/lambda_bc_EJ_server.RData'))
load(here::here('09_simulation_and_approximation-cdf/final_models.RData'))

final_models <- tibble(
    test = rep(c('all', 'e_j'), times = c(3, 3)),
    case = rep(1:3, times = 2),
    models = list(clean_mod_all_1, clean_mod_all_2, clean_mod_all_3,
                  clean_mod_E_J_1, clean_mod_E_J_2, clean_mod_E_J_3)
)

# Load Simulation Data 
if(Sys.info()['nodename'] == "DELL-ARBEIT") { # Jens 
    Data <- readRDS('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_100k.rds')
    # load('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_1_m.RData')
} else if(Sys.info()['user'] == "Janine") { # Janine
    load("/Users/Janine/Desktop/BayerHanck/Data_1_m.RData")
} else if(Sys.info()['nodename'] == "OEK-TS01") { # Server
    load('D:\\Klenke\\Data_1_m.RData')
}

lambda_p <- get_lambda(lambda_bc_ALL, 1, 'p')

# Cut Data 
Data %<>%
    dplyr::select(-p_value_Fisher_all) %>%
    dplyr::mutate(p_value_Fisher_bc = ((p_value_Fisher_E_J^lambda_p)-1)/lambda_p, 
                  p_value_Fisher_lg = log(p_value_Fisher_E_J)) %>%
    dplyr::rename(p_value_Fisher = p_value_Fisher_E_J) %>%
    dplyr::filter(p_value_Fisher %in% c(min(Data$p_value_Fisher_E_J), seq(0, 1, 0.0001))) 

### 


# plot data
data_case_1 <- Data %>%
    dplyr::filter(case == 1) %>%
    dplyr::mutate(stat_Fisher_E_J_bc = ((stat_Fisher_E_J^get_lambda(lambda_bc_EJ, 1, 'stat'))-1)/get_lambda(lambda_bc_EJ, 1, 'stat'),
                  stat_Fisher_all_bc = ((stat_Fisher_all^get_lambda(lambda_bc_ALL, 1, 'stat'))-1)/get_lambda(lambda_bc_ALL, 1, 'stat')) %>%
    modelr::add_predictions(clean_mod_all_1) %>%
    dplyr::mutate(PRED_all = invBoxCox(pred)) %>%
    dplyr::mutate(PRED_cor_all = case_when(
        PRED_all <= 0 ~ 1e-12,
        PRED_all >= 1 ~ 1 - 1e-12,
        TRUE ~ PRED_all))


#-- Analysis  tables ----

sum_table_all_case_1 <- best_5_table(table_all_case_1, 'calls_all')

sum_table_all_case_2 <- best_5_table(table_all_case_2, 'calls_all')

sum_table_all_case_3 <- best_5_table(table_all_case_3, 'calls_all')

sum_table_E_J_case_1 <- best_5_table(table_E_J_case_1, 'calls_E_J')

sum_table_E_J_case_2 <- best_5_table(table_E_J_case_2, 'calls_E_J')

sum_table_E_J_case_3 <- best_5_table(table_E_J_case_3, 'calls_E_J')

#-- Analysis  plots ----

data_case_1 %>%
    dplyr::filter(p_value_Fisher <= 0.2) %>%
    own_plot(max_graph = 0.2)


own_plot(data_case_1)

