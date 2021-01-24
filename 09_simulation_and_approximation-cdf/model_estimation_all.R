#---- Functions ----
# metric function
metric_fun <- function(object){
    
    # save call and case
    mod_call <- paste(object$call$form[2], object$call$form[1], object$call$form[3])
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
        dependent = rep(seq(0+(1/(nrow(object$model)/11)), 1, 1/(nrow(object$model)/11)), 11)
        #object$model$p_value_Fisher_E_J[1:5]
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
                      RMSE = as.numeric(RMSE),
                      RMSE_cor = as.numeric(RMSE_cor),
                      RMSE_0.2 = as.numeric(RMSE_0.2), 
                      RMSE_cor_0.2 = as.numeric(RMSE_cor_0.2), 
                      case = as.character(case),
                      pred = list(values)
    )
    
    return(mod_sum)
}

# bind functionen
bind_model_metrics <- function(new_metrics, old_metrics =  model_metrics_ALL) {
    model_metrics_all <- rbind(old_metrics,
                               new_metrics)
    
    return(model_metrics_all)
}

# inverse BoxCox function
invBoxCox <- function(x){
    x <- if (lambda_p == 0) exp(as.complex(x)) else (lambda_p*as.complex(x) + 1)^(1/lambda_p)
    return(Re(x))
}

# delete function 
delete_fun <- function(){
    rm(list=setdiff(ls(.GlobalEnv), important_things), envir = .GlobalEnv)
}

#---- Preliminary ---- 
source(here::here('01_code/packages/packages.R'))

## Load Simulation Data 
if(Sys.info()['nodename'] == "DESKTOP-ATT92OH") { # Jens 
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

#---- Boxcox Transformation ----

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
    case = rep(1:3, each = 2), 
    side = rep(c('stat', 'p'), times = 3),
    value = c(lambda_stat_case_1, lambda_p, lambda_stat_case_2, lambda_p,  lambda_stat_case_3, lambda_p)   
)

#---- Set up metrics -----

model_metrics_ALL <- NULL

# save important things
important_things <- c(ls(), 
                      "important_things")


# ----  Models - Case 1 ----
tictoc::tic()

models_poly <- list(
    # functional form: poly(t, p) + (1/k)
    mod_all_case.1_p_3 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + I(1/k),
                             data = data_case_1),
    
    mod_all_case.1_p_4 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + I(1/k),
                             data = data_case_1),
    
    mod_all_case.1_p_5 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + I(1/k),
                             data = data_case_1),
    
    mod_all_case.1_p_6 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + I(1/k),
                             data = data_case_1), 
    
    # functional form: poly(t, p) + (1/k) + poly(t, p)*(1/k)
    
    mod_all_case.1_p_3_3 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + I(1/k) + poly(k, 3)*(1/k),
                               data = data_case_1),
    
    mod_all_case.1_p_4_4 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + I(1/k) + poly(k, 4)*(1/k),
                               data = data_case_1),
    
    mod_all_case.1_p_5_5 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + I(1/k) + poly(k, 5)*(1/k),
                               data = data_case_1),
    
    mod_all_case.1_p_6_6 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + I(1/k) + poly(k, 6)*(1/k),
                               data = data_case_1)
)

# save metrics
model_metrics_ALL <- bind_model_metrics(models_poly %>% purrr::map_dfr(metric_fun))


# delete model
delete_fun()

models_poly_2 <- list(
    # functional form: poly(t, p) + log(k) + poly(t, p)*log(k)
    mod_all_case.1_p_3_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + log(k) + poly(k, 3)*log(k),
                                 data = data_case_1),
    
    mod_all_case.1_p_4_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + log(k) + poly(k, 4)*log(k),
                                 data = data_case_1),
    
    mod_all_case.1_p_5_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + log(k) + poly(k, 5)*log(k),
                                 data = data_case_1),
    
    mod_all_case.1_p_6_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + log(k) + poly(k, 6)*log(k),
                                 data = data_case_1),
    
    # functional form: poly(t, p) + k + (1/k)
    mod_all_case.1_p_3_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + k + I(1/k),
                                 data = data_case_1),
    
    mod_all_case.1_p_4_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + k + I(1/k),
                                 data = data_case_1),
    
    mod_all_case.1_p_5_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + k + I(1/k),
                                 data = data_case_1),
    
    mod_all_case.1_p_6_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + k + I(1/k),
                                 data = data_case_1)
)

# save metrics
model_metrics_ALL <- bind_model_metrics(models_poly_2 %>% purrr::map_dfr(metric_fun))

# delete model
delete_fun()

# log()
models_log <- list(
    mod_all_case.1_log_1 = lm(p_value_Fisher_all ~ log(stat_Fisher_all) + k,
                              data = data_case_1),
    
    mod_all_case.1_log_2 = lm(p_value_Fisher_all ~ log(stat_Fisher_all) * k,
                              data = data_case_1),
    
    mod_all_case.1_poly_log_3 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 3) + log(k),
                                    data = data_case_1),
    
    mod_all_case.1_poly_log_4 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 4) + log(k),
                                    data = data_case_1),
    
    mod_all_case.1_poly_log_5 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 5) + log(k),
                                    data = data_case_1),
    
    mod_all_case.1_poly_log_6 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 6) + log(k),
                                    data = data_case_1),
    
    mod_all_case.1_poly_log_7 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 7) + log(k),
                                    data = data_case_1),
    
    mod_all_case.1_poly_log_8 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 8) + log(k),
                                    data = data_case_1), 
    
    mod_all_case.1_poly_log_9 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 9) + log(k),
                                    data = data_case_1),
    
    mod_all_case.1_poly_log_10 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 10) + log(k),
                                     data = data_case_1),
    
    mod_all_case.1_poly_log_m_3 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 3) * log(k),
                                      data = data_case_1),
    
    mod_all_case.1_poly_log_m_4 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 4) * log(k),
                                      data = data_case_1),
    
    mod_all_case.1_poly_log_m_5 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 5) * log(k),
                                      data = data_case_1),
    
    mod_all_case.1_poly_log_m_6 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 6)  *log(k),
                                      data = data_case_1),
    
    mod_all_case.1_poly_log_m_7 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 7) * log(k),
                                      data = data_case_1),
    
    mod_all_case.1_poly_log_m_8 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 8) * log(k),
                                      data = data_case_1), 
    
    mod_all_case.1_poly_log_m_9 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 9) * log(k),
                                      data = data_case_1),
    
    mod_all_case.1_poly_log_m_10 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 10) * log(k),
                                       data = data_case_1)
)

# save metrics
model_metrics_ALL <- bind_model_metrics(models_log %>% purrr::map_dfr(metric_fun))

# delete model
delete_fun()


models_bc <- list(
    mod_all_case.1_bc_poly_log_m_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k),
                                          data = data_case_1),
    
    mod_all_case.1_bc_poly_m_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * k, 
                                      data = data_case_1),
    
    mod_all_case.1_bc_poly_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) + k, 
                                    data = data_case_1),
    
    mod_all_case.1_bc_poly_log_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) + log(k), 
                                        data = data_case_1),
    
    mod_all_case.1_bc_poly_log_m_10_I <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + I(1/k),
                                            data = data_case_1),
    
    mod_all_case.1_bc_poly_m_10_I <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * k + I(1/k),
                                        data = data_case_1),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k),
                                                   data = data_case_1),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_sqrt <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k) + sqrt(k),
                                                        data = data_case_1),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_sqrt <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*sqrt(k),
                                                      data = data_case_1)
)

model_metrics_ALL <- bind_model_metrics(models_bc %>% purrr::map_dfr(metric_fun))

delete_fun()

# bc y
models_bc_y <- list(
    mod_all_case.1_bc_poly_log_m_6_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6)*log(k), 
                                           data = data_case_1),
    
    mod_all_case.1_bc_poly_log_m_10_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10)*log(k), 
                                            data = data_case_1), 
    
    mod_all_case.1_bc_poly_log_6_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6) + log(k), 
                                         data = data_case_1), 
    
    mod_all_case.1_bc_poly_log_10_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10) + log(k), 
                                          data = data_case_1), 
    
    mod_all_case.1_bc_poly_log_m_6_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6)*log(k) + I(1/k), 
                                             data = data_case_1), 
    
    mod_all_case.1_bc_poly_log_m_10_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10)*log(k) + I(1/k), 
                                              data = data_case_1),
    
    mod_all_case.1_bc_poly_log_m_6_poly_m_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6)*log(k) + poly(stat_Fisher_all_bc, 6)*I(1/k), 
                                                    data = data_case_1),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10)*log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k), 
                                                     data = data_case_1), 
    
    mod_all_case.1_bc_bc <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*sqrt(k),
                               data = data_case_1)
)

model_metrics_ALL <- bind_model_metrics(models_bc_y %>% purrr::map_dfr(metric_fun))

delete_fun()

# log y
models_exp <- list(
    mod_all_case.1_bc_poly_log_m_10_I_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10)*log(k) + I(1/k), 
                                               data = data_case_1),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10)*log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k), 
                                                      data = data_case_1),
    
    mod_all_case.1_bc_poly_log_m_10_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10)*log(k), 
                                             data = data_case_1), 
    
    mod_all_case.1_bc_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*sqrt(k),
                               data = data_case_1),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_sqrt_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k) + sqrt(k),
                                                           data = data_case_1)
)

model_metrics_ALL <- bind_model_metrics(models_exp %>% purrr::map_dfr(metric_fun))

delete_fun()

tictoc::toc()

# ----  Models - Case 2 ----

tictoc::tic('case_2')

models_poly <- list(
    # functional form: poly(t, p) + (1/k)
    mod_all_case.1_p_3 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + I(1/k),
                             data = data_case_2),
    
    mod_all_case.1_p_4 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + I(1/k),
                             data = data_case_2),
    
    mod_all_case.1_p_5 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + I(1/k),
                             data = data_case_2),
    
    mod_all_case.1_p_6 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + I(1/k),
                             data = data_case_2), 
    
    # functional form: poly(t, p) + (1/k) + poly(t, p)*(1/k)
    
    mod_all_case.1_p_3_3 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + I(1/k) + poly(k, 3)*(1/k),
                               data = data_case_2),
    
    mod_all_case.1_p_4_4 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + I(1/k) + poly(k, 4)*(1/k),
                               data = data_case_2),
    
    mod_all_case.1_p_5_5 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + I(1/k) + poly(k, 5)*(1/k),
                               data = data_case_2),
    
    mod_all_case.1_p_6_6 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + I(1/k) + poly(k, 6)*(1/k),
                               data = data_case_2)
)

# save metrics
model_metrics_ALL <- bind_model_metrics(models_poly %>% purrr::map_dfr(metric_fun))


# delete model
delete_fun()

models_poly_2 <- list(
    # functional form: poly(t, p) + log(k) + poly(t, p)*log(k)
    mod_all_case.1_p_3_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + log(k) + poly(k, 3)*log(k),
                                 data = data_case_2),
    
    mod_all_case.1_p_4_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + log(k) + poly(k, 4)*log(k),
                                 data = data_case_2),
    
    mod_all_case.1_p_5_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + log(k) + poly(k, 5)*log(k),
                                 data = data_case_2),
    
    mod_all_case.1_p_6_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + log(k) + poly(k, 6)*log(k),
                                 data = data_case_2),
    
    # functional form: poly(t, p) + k + (1/k)
    mod_all_case.1_p_3_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + k + I(1/k),
                                 data = data_case_2),
    
    mod_all_case.1_p_4_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + k + I(1/k),
                                 data = data_case_2),
    
    mod_all_case.1_p_5_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + k + I(1/k),
                                 data = data_case_2),
    
    mod_all_case.1_p_6_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + k + I(1/k),
                                 data = data_case_2)
)

# save metrics
model_metrics_ALL <- bind_model_metrics(models_poly_2 %>% purrr::map_dfr(metric_fun))

# delete model
delete_fun()

# log()
models_log <- list(
    mod_all_case.1_log_1 = lm(p_value_Fisher_all ~ log(stat_Fisher_all) + k,
                              data = data_case_2),
    
    mod_all_case.1_log_2 = lm(p_value_Fisher_all ~ log(stat_Fisher_all) * k,
                              data = data_case_2),
    
    mod_all_case.1_poly_log_3 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 3) + log(k),
                                    data = data_case_2),
    
    mod_all_case.1_poly_log_4 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 4) + log(k),
                                    data = data_case_2),
    
    mod_all_case.1_poly_log_5 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 5) + log(k),
                                    data = data_case_2),
    
    mod_all_case.1_poly_log_6 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 6) + log(k),
                                    data = data_case_2),
    
    mod_all_case.1_poly_log_7 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 7) + log(k),
                                    data = data_case_2),
    
    mod_all_case.1_poly_log_8 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 8) + log(k),
                                    data = data_case_2), 
    
    mod_all_case.1_poly_log_9 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 9) + log(k),
                                    data = data_case_2),
    
    mod_all_case.1_poly_log_10 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 10) + log(k),
                                     data = data_case_2),
    
    mod_all_case.1_poly_log_m_3 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 3) * log(k),
                                      data = data_case_2),
    
    mod_all_case.1_poly_log_m_4 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 4) * log(k),
                                      data = data_case_2),
    
    mod_all_case.1_poly_log_m_5 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 5) * log(k),
                                      data = data_case_2),
    
    mod_all_case.1_poly_log_m_6 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 6)  *log(k),
                                      data = data_case_2),
    
    mod_all_case.1_poly_log_m_7 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 7) * log(k),
                                      data = data_case_2),
    
    mod_all_case.1_poly_log_m_8 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 8) * log(k),
                                      data = data_case_2), 
    
    mod_all_case.1_poly_log_m_9 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 9) * log(k),
                                      data = data_case_2),
    
    mod_all_case.1_poly_log_m_10 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 10) * log(k),
                                       data = data_case_2)
)

# save metrics
model_metrics_ALL <- bind_model_metrics(models_log %>% purrr::map_dfr(metric_fun))

# delete model
delete_fun()

models_bc <- list(
    mod_all_case.1_bc_poly_log_m_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k),
                                          data = data_case_2),
    
    mod_all_case.1_bc_poly_m_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * k, 
                                      data = data_case_2),
    
    mod_all_case.1_bc_poly_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) + k, 
                                    data = data_case_2),
    
    mod_all_case.1_bc_poly_log_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) + log(k), 
                                        data = data_case_2),
    
    mod_all_case.1_bc_poly_log_m_10_I <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + I(1/k),
                                            data = data_case_2),
    
    mod_all_case.1_bc_poly_m_10_I <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * k + I(1/k),
                                        data = data_case_2),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k),
                                                   data = data_case_2),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_sqrt <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k) + sqrt(k),
                                                        data = data_case_2),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_sqrt <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*sqrt(k),
                                                      data = data_case_2)
)

model_metrics_ALL <- bind_model_metrics(models_bc %>% purrr::map_dfr(metric_fun))

delete_fun()

# bc y
models_bc_y <- list(
    mod_all_case.1_bc_poly_log_m_6_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6)*log(k), 
                                           data = data_case_2), 
    
    mod_all_case.1_bc_poly_log_m_10_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10)*log(k), 
                                            data = data_case_2), 
    
    mod_all_case.1_bc_poly_log_6_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6) + log(k), 
                                         data = data_case_2), 
    
    mod_all_case.1_bc_poly_log_10_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10) + log(k), 
                                          data = data_case_2), 
    
    mod_all_case.1_bc_poly_log_m_6_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6)*log(k) + I(1/k), 
                                             data = data_case_2), 
    
    mod_all_case.1_bc_poly_log_m_10_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10)*log(k) + I(1/k), 
                                              data = data_case_2),
    
    mod_all_case.1_bc_poly_log_m_6_poly_m_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6)*log(k) + poly(stat_Fisher_all_bc, 6)*I(1/k), 
                                                    data = data_case_2),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10)*log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k), 
                                                     data = data_case_2), 
    
    mod_all_case.1_bc_bc <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*sqrt(k),
                               data = data_case_2)
)

model_metrics_ALL <- bind_model_metrics(models_bc_y %>% purrr::map_dfr(metric_fun))

delete_fun()

# log y
models_exp <- list(
    mod_all_case.1_bc_poly_log_m_10_I_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10)*log(k) + I(1/k), 
                                               data = data_case_2),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10)*log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k), 
                                                      data = data_case_2),
    
    mod_all_case.1_bc_poly_log_m_10_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10)*log(k), 
                                             data = data_case_2), 
    
    mod_all_case.1_bc_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*sqrt(k),
                               data = data_case_2),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_sqrt_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k) + sqrt(k),
                                                           data = data_case_2)
)

model_metrics_ALL <- bind_model_metrics(models_exp %>% purrr::map_dfr(metric_fun))

delete_fun()

tictoc::toc()

# ----  Models - Case 3 ----

tictoc::tic()

models_poly <- list(
    # functional form: poly(t, p) + (1/k)
    mod_all_case.1_p_3 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + I(1/k),
                             data = data_case_3),
    
    mod_all_case.1_p_4 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + I(1/k),
                             data = data_case_3),
    
    mod_all_case.1_p_5 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + I(1/k),
                             data = data_case_3),
    
    mod_all_case.1_p_6 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + I(1/k),
                             data = data_case_3), 
    
    # functional form: poly(t, p) + (1/k) + poly(t, p)*(1/k)
    
    mod_all_case.1_p_3_3 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + I(1/k) + poly(k, 3)*(1/k),
                               data = data_case_3),
    
    mod_all_case.1_p_4_4 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + I(1/k) + poly(k, 4)*(1/k),
                               data = data_case_3),
    
    mod_all_case.1_p_5_5 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + I(1/k) + poly(k, 5)*(1/k),
                               data = data_case_3),
    
    mod_all_case.1_p_6_6 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + I(1/k) + poly(k, 6)*(1/k),
                               data = data_case_3)
)

# save metrics
model_metrics_ALL <- bind_model_metrics(models_poly %>% purrr::map_dfr(metric_fun))


# delete model
delete_fun()

models_poly_2 <- list(
    # functional form: poly(t, p) + log(k) + poly(t, p)*log(k)
    mod_all_case.1_p_3_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + log(k) + poly(k, 3)*log(k),
                                 data = data_case_3),
    
    mod_all_case.1_p_4_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + log(k) + poly(k, 4)*log(k),
                                 data = data_case_3),
    
    mod_all_case.1_p_5_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + log(k) + poly(k, 5)*log(k),
                                 data = data_case_3),
    
    mod_all_case.1_p_6_log <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + log(k) + poly(k, 6)*log(k),
                                 data = data_case_3),
    
    # functional form: poly(t, p) + k + (1/k)
    mod_all_case.1_p_3_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 3) + k + I(1/k),
                                 data = data_case_3),
    
    mod_all_case.1_p_4_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 4) + k + I(1/k),
                                 data = data_case_3),
    
    mod_all_case.1_p_5_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 5) + k + I(1/k),
                                 data = data_case_3),
    
    mod_all_case.1_p_6_k_1 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all, 6) + k + I(1/k),
                                 data = data_case_3)
)

# save metrics
model_metrics_ALL <- bind_model_metrics(models_poly_2 %>% purrr::map_dfr(metric_fun))

# delete model
delete_fun()

# log()
models_log <- list(
    mod_all_case.1_log_1 = lm(p_value_Fisher_all ~ log(stat_Fisher_all) + k,
                              data = data_case_3),
    
    mod_all_case.1_log_2 = lm(p_value_Fisher_all ~ log(stat_Fisher_all) * k,
                              data = data_case_3),
    
    mod_all_case.1_poly_log_3 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 3) + log(k),
                                    data = data_case_3),
    
    mod_all_case.1_poly_log_4 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 4) + log(k),
                                    data = data_case_3),
    
    mod_all_case.1_poly_log_5 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 5) + log(k),
                                    data = data_case_3),
    
    mod_all_case.1_poly_log_6 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 6) + log(k),
                                    data = data_case_3),
    
    mod_all_case.1_poly_log_7 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 7) + log(k),
                                    data = data_case_3),
    
    mod_all_case.1_poly_log_8 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 8) + log(k),
                                    data = data_case_3), 
    
    mod_all_case.1_poly_log_9 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 9) + log(k),
                                    data = data_case_3),
    
    mod_all_case.1_poly_log_10 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 10) + log(k),
                                     data = data_case_3),
    
    mod_all_case.1_poly_log_m_3 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 3) * log(k),
                                      data = data_case_3),
    
    mod_all_case.1_poly_log_m_4 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 4) * log(k),
                                      data = data_case_3),
    
    mod_all_case.1_poly_log_m_5 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 5) * log(k),
                                      data = data_case_3),
    
    mod_all_case.1_poly_log_m_6 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 6)  *log(k),
                                      data = data_case_3),
    
    mod_all_case.1_poly_log_m_7 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 7) * log(k),
                                      data = data_case_3),
    
    mod_all_case.1_poly_log_m_8 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 8) * log(k),
                                      data = data_case_3), 
    
    mod_all_case.1_poly_log_m_9 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 9) * log(k),
                                      data = data_case_3),
    
    mod_all_case.1_poly_log_m_10 <- lm(p_value_Fisher_all ~ poly(log(stat_Fisher_all), 10) * log(k),
                                       data = data_case_3)
)

# save metrics
model_metrics_ALL <- bind_model_metrics(models_log %>% purrr::map_dfr(metric_fun))

# delete model
delete_fun()


models_bc <- list(
    mod_all_case.1_bc_poly_log_m_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k),
                                          data = data_case_3),
    
    mod_all_case.1_bc_poly_m_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * k, 
                                      data = data_case_3),
    
    mod_all_case.1_bc_poly_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) + k, 
                                    data = data_case_3),
    
    mod_all_case.1_bc_poly_log_10 <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) + log(k), 
                                        data = data_case_3),
    
    mod_all_case.1_bc_poly_log_m_10_I <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + I(1/k),
                                            data = data_case_3),
    
    mod_all_case.1_bc_poly_m_10_I <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * k + I(1/k),
                                        data = data_case_3),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k),
                                                   data = data_case_3),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_sqrt <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k) + sqrt(k),
                                                        data = data_case_3),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_sqrt <- lm(p_value_Fisher_all ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*sqrt(k),
                                                      data = data_case_3)
)

model_metrics_ALL <- bind_model_metrics(models_bc %>% purrr::map_dfr(metric_fun))

delete_fun()


# bc y
models_bc_y <- list(
    mod_all_case.1_bc_poly_log_m_6_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6)*log(k), 
                                           data = data_case_3), 
    
    mod_all_case.1_bc_poly_log_m_10_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10)*log(k), 
                                            data = data_case_3), 
    
    mod_all_case.1_bc_poly_log_6_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6) + log(k), 
                                         data = data_case_3), 
    
    mod_all_case.1_bc_poly_log_10_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10) + log(k), 
                                          data = data_case_3), 
    
    mod_all_case.1_bc_poly_log_m_6_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6)*log(k) + I(1/k), 
                                             data = data_case_3), 
    
    mod_all_case.1_bc_poly_log_m_10_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10)*log(k) + I(1/k), 
                                              data = data_case_3),
    
    mod_all_case.1_bc_poly_log_m_6_poly_m_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 6)*log(k) + poly(stat_Fisher_all_bc, 6)*I(1/k), 
                                                    data = data_case_3),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_y <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10)*log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k), 
                                                     data = data_case_3), 
    
    mod_all_case.1_bc_bc <- lm(p_value_Fisher_all_bc ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*sqrt(k),
                               data = data_case_3)
)

model_metrics_ALL <- bind_model_metrics(models_bc_y %>% purrr::map_dfr(metric_fun))

delete_fun()

# log y
models_exp <- list(
    mod_all_case.1_bc_poly_log_m_10_I_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10)*log(k) + I(1/k), 
                                               data = data_case_3),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10)*log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k), 
                                                      data = data_case_3),
    
    mod_all_case.1_bc_poly_log_m_10_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10)*log(k), 
                                             data = data_case_3), 
    
    mod_all_case.1_bc_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*sqrt(k),
                               data = data_case_3),
    
    mod_all_case.1_bc_poly_log_m_10_poly_m_I_sqrt_lg <- lm(p_value_Fisher_all_lg ~ poly(stat_Fisher_all_bc, 10) * log(k) + poly(stat_Fisher_all_bc, 10)*I(1/k) + sqrt(k),
                                                           data = data_case_3)
)

model_metrics_ALL <- bind_model_metrics(models_exp %>% purrr::map_dfr(metric_fun))

delete_fun()

tictoc::toc()

#---- saving metrics ----

save(model_metrics_ALL, file = here::here("09_simulation_and_approximation-cdf/model_metrics_ALL.Rdata"))
save(lambda_bc_ALL, file = here::here("09_simulation_and_approximation-cdf/lambda_bc_ALL.Rdata"))
