#-- functions ----
clean_model <- function(object){
    object$fitted.values <- NULL
    object$residuals <- NULL
    object$effects <- NULL
    object$model <- NULL
    object$qr$qr <- NULL 
    
    return(object)
}

#---- Preliminary ---- 
source(here::here('01_code/packages/packages.R'))

# Load Simulation Data 

# Data
if(Sys.info()['nodename'] == "DELL-ARBEIT") { # Jens 
    Data <- readRDS('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_100k.rds')
    # load('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_1_m.RData')
} else if(Sys.info()['nodename'] == "MacBook-Pro.local") { # Janine
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

# case_1 
lambda_stat_case_1 <- Rfast::bc(data_case_1$stat_Fisher_all)
lambda_p <- Rfast::bc(data_case_1$p_value_Fisher_all)

data_case_1 <- data_case_1 %>% mutate(
    stat_Fisher_all_bc = ((stat_Fisher_all^lambda_stat_case_1)-1)/lambda_stat_case_1,
    p_value_Fisher_all_bc = ((p_value_Fisher_all^lambda_p)-1)/lambda_p, 
    p_value_Fisher_all_lg = log(p_value_Fisher_all)
)

# case_2 
lambda_stat_case_2 <- Rfast::bc(data_case_2$stat_Fisher_all)

data_case_2 <- data_case_2 %>% mutate(
    stat_Fisher_all_bc = ((stat_Fisher_all^lambda_stat_case_2)-1)/lambda_stat_case_2,
    p_value_Fisher_all_bc = ((p_value_Fisher_all^lambda_p)-1)/lambda_p,
    p_value_Fisher_all_lg = log(p_value_Fisher_all)
)

# case_3 
lambda_stat_case_3 <- Rfast::bc(data_case_3$stat_Fisher_all)

data_case_3 <- data_case_3 %>% mutate(
    stat_Fisher_all_bc = ((stat_Fisher_all^lambda_stat_case_3)-1)/lambda_stat_case_3,
    p_value_Fisher_all_bc = ((p_value_Fisher_all^lambda_p)-1)/lambda_p,
    p_value_Fisher_all_lg = log(p_value_Fisher_all)
)

lambda_bc_ALL <- tibble(
    test = rep("all", 3),
    case = 1:3, 
    bc_p = c(lambda_p, lambda_p, lambda_p),
    bc_stat = c(lambda_stat_case_1, lambda_stat_case_2,  lambda_stat_case_3)   
)

# case_1 
lambda_stat_case_1 <- Rfast::bc(data_case_1$stat_Fisher_E_J)
lambda_p <- Rfast::bc(data_case_1$p_value_Fisher_E_J)

data_case_1 <- data_case_1 %>% mutate(
    stat_Fisher_E_J_bc = ((stat_Fisher_E_J^lambda_stat_case_1)-1)/lambda_stat_case_1,
    p_value_Fisher_E_J_bc = ((p_value_Fisher_E_J^lambda_p)-1)/lambda_p, 
    p_value_Fisher_E_J_lg = log(p_value_Fisher_E_J)
)

# case_2 
lambda_stat_case_2 <- Rfast::bc(data_case_2$stat_Fisher_E_J)

data_case_2 <- data_case_2 %>% mutate(
    stat_Fisher_E_J_bc = ((stat_Fisher_E_J^lambda_stat_case_2)-1)/lambda_stat_case_2,
    p_value_Fisher_E_J_bc = ((p_value_Fisher_E_J^lambda_p)-1)/lambda_p,
    p_value_Fisher_E_J_lg = log(p_value_Fisher_E_J)
)

# case_3 
lambda_stat_case_3 <- Rfast::bc(data_case_3$stat_Fisher_E_J)

data_case_3 <- data_case_3 %>% mutate(
    stat_Fisher_E_J_bc = ((stat_Fisher_E_J^lambda_stat_case_3)-1)/lambda_stat_case_3,
    p_value_Fisher_E_J_bc = ((p_value_Fisher_E_J^lambda_p)-1)/lambda_p,
    p_value_Fisher_E_J_lg = log(p_value_Fisher_E_J)
)

lambda_bc_EJ <-  tibble(
        test = rep("eg-j", 3),
        case = 1:3, 
        bc_p = c(lambda_p, lambda_p, lambda_p),
        bc_stat = c(lambda_stat_case_1, lambda_stat_case_2,  lambda_stat_case_3)   
    )
    
#-- set up final bc ----

final_model_metrics <- bind_rows(lambda_bc_ALL, lambda_bc_EJ)


#metrics
load(here::here("09_simulation_and_approximation-cdf/model_metrics_ALL_server.Rdata"))

load(here::here("09_simulation_and_approximation-cdf/model_metrics_E_J_server.Rdata"))

#-- Data models ---- 
best_models <- bind_rows( 
    model_metrics_E_J%>%
        dplyr::group_by(case)%>%
        dplyr::filter(RMSE_cor_0.2 == min(RMSE_cor_0.2)) # trade off maybe we should use RMSE_0.2_cor
    ,
    
    model_metrics_ALL%>%
        dplyr::group_by(case)%>%
        dplyr::filter(RMSE_cor_0.2 == min(RMSE_cor_0.2)) # trade off maybe we should use RMSE_0.2_cor
)%>%
    dplyr::select(model, case)


#-- Fitting Final models ----

# E_J
mod_E_J_case_1 <- clean_model(lm(as.formula(best_models$model[1]) , data = data_case_1))

mod_E_J_case_2 <- clean_model(lm(as.formula(best_models$model[2]) , data = data_case_2))

mod_E_J_case_3 <- clean_model(lm(as.formula(best_models$model[3]) , data = data_case_3))

# ALL
mod_ALL_case_1 <- clean_model(lm(as.formula(best_models$model[4]) , data = data_case_1))

mod_ALL_case_2 <- clean_model(lm(as.formula(best_models$model[5]) , data = data_case_2))

mod_ALL_case_3 <- clean_model(lm(as.formula(best_models$model[6]) , data = data_case_3))

#-- save final models ---
save(final_model_metrics, mod_E_J_case_1, mod_E_J_case_2, mod_E_J_case_3, mod_ALL_case_1,
     mod_ALL_case_2, mod_ALL_case_3, 
     file = here::here("09_simulation_and_approximation-cdf/p_value_approx.RData"))
