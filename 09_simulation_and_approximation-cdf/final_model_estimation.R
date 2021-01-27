#-- functions ----
clean_model <- function(object){
    object$fitted.values <- NULL
    object$residuals <- NULL
    object$effects <- NULL
    object$model <- NULL
    object$qr$qr <- NULL 
    
    return(object)
}

# metric function
metric_fun <- function(object){
    
    # save call and case
    mod_call <- deparse(substitute(object))
    case <- object$call$data
    
    # dependent variable 
    dep_var <- as.character(object$call$form[2])
    
    # dataset for prediction
    values <- tibble(
        PRED = if(str_sub(dep_var , nchar(dep_var)-2, nchar(dep_var)) == '_bc'){ 
            invBoxCox(object$fitted.values)
        } else if(str_sub(dep_var , nchar(dep_var)-2, nchar(dep_var)) == '_lg'){
            exp(object$fitted.values)
        } else {object$fitted.values}, # y_hat 
        dependent = rep(seq(0+(1/(nrow(object$model)/11)), 1, 1/(nrow(object$model)/11)), 11),
        k = rep(1:11, each = nrow(object$model)/11)
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
        dplyr::filter(dependent >= 0.8)
    RMSE_0.2 <- sqrt((sum((values_0.2$PRED - values_0.2$dependent)^2)/nrow(values_0.2)))
    RMSE_cor_0.2 <- sqrt((sum((values_0.2$PRED_cor - values_0.2$dependent)^2)/nrow(values_0.2)))
    
    # values for the plots 
    values <- values%>%
        dplyr::filter(dependent %in% seq(0.001, 1, 0.001))
    
    mod_sum <- tibble(model = as.character(mod_call),
                      case = as.character(case),
                      pred = list(values)
    )
    
    return(mod_sum)
}

# bind functionen
bind_model_metrics <- function(new_metrics, old_metrics =  model_metrics) {
    model_metrics <- rbind(old_metrics,
                               new_metrics)
    
    return(model_metrics)
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
    model_metrics_ALL%>%
        dplyr::filter(case == 'data_case_1',
                      model == 'p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10) * sqrt(k)'),
    
    model_metrics_ALL%>%
        dplyr::filter(case == 'data_case_2',
                      model == 'p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10) * sqrt(k)'),
    
    model_metrics_ALL%>%
        dplyr::filter(case == 'data_case_3',
                      model == 'p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10) * sqrt(k)'),
    
    model_metrics_E_J%>%
        dplyr::filter(case == 'data_case_1',
                        model == 'p_value_Fisher_E_J_lg ~ poly(stat_Fisher_E_J_bc, 10) * log(k) + poly(stat_Fisher_E_J_bc, 10) * I(1/k)'),
    
    model_metrics_E_J%>%
        dplyr::filter(case == 'data_case_2',
                        model == 'p_value_Fisher_E_J_lg ~ poly(stat_Fisher_E_J_bc, 10) * log(k)'),
    
    model_metrics_E_J%>%
        dplyr::filter(case == 'data_case_3',
                        model == 'p_value_Fisher_E_J_lg ~ poly(stat_Fisher_E_J_bc, 10) * log(k) + poly(stat_Fisher_E_J_bc, 10) * sqrt(k)'))%>%
    dplyr::select(model, case)


model_metrics <- NULL
#-- Fitting Final models ----

# E_J
mod_E_J_case_1 <- lm(as.formula(best_models$model[1]) , data = data_case_1)

model_metrics <- bind_model_metrics(metric_fun(mod_E_J_case_1))

mod_E_J_case_1 <- clean_model(mod_E_J_case_1)

mod_E_J_case_2 <- lm(as.formula(best_models$model[2]) , data = data_case_2)

model_metrics <- bind_model_metrics(metric_fun(mod_E_J_case_2))

mod_E_J_case_2 <- clean_model(mod_E_J_case_2)

mod_E_J_case_3 <- lm(as.formula(best_models$model[3]) , data = data_case_3)

model_metrics <- bind_model_metrics(metric_fun(mod_E_J_case_3))

mod_E_J_case_3 <- clean_model(mod_E_J_case_3)


# ALL
mod_ALL_case_1 <- lm(as.formula(best_models$model[4]) , data = data_case_1)

model_metrics <- bind_model_metrics(metric_fun(mod_ALL_case_1))

mod_ALL_case_1 <- clean_model(mod_ALL_case_1)

mod_ALL_case_2 <- clean_model(lm(as.formula(best_models$model[5]) , data = data_case_2))

model_metrics <- bind_model_metrics(metric_fun(mod_ALL_case_2))

mod_ALL_case_2 <- clean_model(mod_ALL_case_2)

mod_ALL_case_3 <- clean_model(lm(as.formula(best_models$model[6]) , data = data_case_3))

model_metrics <- bind_model_metrics(metric_fun(mod_ALL_case_3))

mod_ALL_case_3 <- clean_model(mod_ALL_case_3)

model_metrics <- model_metrics%>%
    dplyr::mutate(test = rep(c('eg-j', 'all'), each = 3))

#-- save final models ---
save(final_model_metrics, mod_E_J_case_1, mod_E_J_case_2, mod_E_J_case_3, mod_ALL_case_1,
     mod_ALL_case_2, mod_ALL_case_3, 
     file = here::here("09_simulation_and_approximation-cdf/p_value_approx.RData"))
