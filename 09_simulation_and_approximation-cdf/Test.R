#---- Data ----

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

# Split Dataset in Cases
data_case_1 <- Data %>%
    dplyr::filter(case == 1)
data_case_2 <- Data %>%
    dplyr::filter(case == 2)
data_case_3 <- Data %>%
    dplyr::filter(case == 3)

#---- Boxcox Transformation ----

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

lambda_values <- bind_rows(lambda_bc_all, 
                           lambda_bc_EJ) %>% 
    dplyr::mutate(test = c(rep("all", 6), 
                           rep("e_j", 6)))

#---- Model ----

mod_neu <- lm(p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, 10) * log(k) * I(1/k) * k * sqrt(k), 
              data = data_case_1)


own_pred <- function(mod, data){
    # lm fit
    # dependent variable 
    dep_var <- "p_value_Fisher_bc"
    # fitted values
    fitted_values <- mod$fitted.values
    # cleand model
    mod <- clean_lm(mod)
    
    values <- tibble(
        PRED = if(str_detect(dep_var, '_bc')){ 
            invBoxCox(fitted_values)
        } else if(str_detect(dep_var, '_lg')){
            exp(fitted_values)
        } else {fitted_values}, # y_hat 
        dependent = data$p_value_Fisher,
        k = data$k
    )
    
    # computing corrected predictions
    values <- values%>%
        dplyr::mutate(PRED_cor = case_when(
            PRED <= 0 ~ 1e-12,
            PRED >= 1 ~ 1 - 1e-12,
            TRUE ~ PRED))
    
    # RMSE and corrected RMSE on the full dataset
    RMSE <- sqrt((sum((values$PRED - values$dependent)^2)/nrow(values)))
    RMSE_cor <- sqrt((sum((values$PRED_cor - values$dependent)^2)/nrow(values)))
    
    # smaller dataset only the intersting part, correction is at 
    values_0.2 <- values %>%
        dplyr::filter(dependent >= 0.2)
    RMSE_0.2 <- sqrt((sum((values_0.2$PRED - values_0.2$dependent)^2)/nrow(values_0.2)))
    RMSE_cor_0.2 <- sqrt((sum((values_0.2$PRED_cor - values_0.2$dependent)^2)/nrow(values_0.2)))
    
    mod_sum <- tibble(model  = list(mod),
                      RMSE = as.numeric(RMSE),
                      RMSE_cor = as.numeric(RMSE_cor),
                      RMSE_0.2 = as.numeric(RMSE_0.2), 
                      RMSE_cor_0.2 = as.numeric(RMSE_cor_0.2) 
    )
    
    return(mod_sum)
}

mod_neu_clean <- mod_neu
mod_neu_clean$qr <- NULL

object_size(mod_neu_clean)

own_pred(mod_neu, data_case_1)
# RMSE_cor_0.2: 0.000466

# Hier habe ich dann alles aus model_evaluation geladen 

models$models[[1]] <- mod_neu

data_case_1_all <- plot_data(Data, 1, 'all')

### case 1
# all
own_plot(data_case_1_all)

# anderes modell
data_case_1_all %>%
    dplyr::filter(p_value_Fisher <= 0.2) %>%
    own_plot(max_graph = 0.2)

alt_coef <- clean_mod_all_1$coefficients %>% as.numeric()
neu_coef <- mod_neu$coefficients

comp <- tibble(alt = alt_coef, neu = neu_coef)

tictoc::tic()
mod_neu_fast <- RcppEigen::fastLm(p_value_Fisher_bc ~ 1 + poly(stat_Fisher_all_bc, 10) * log(k) * I(1/k) * k * sqrt(k), 
                    data = data_case_1)
tictoc::toc()
# 2148.642 sec elapsed
object_size(mod_neu_fast)
# 176 MB

