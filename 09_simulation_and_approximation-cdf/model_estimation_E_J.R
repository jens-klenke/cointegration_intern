#-----------------------
# Functions
#-----------------------
metric_fun <- function(object){
    
    # dependent variable 
    dep_var <- as.character(object$call$form[2])
    
    # interesting border
    dep_var_int_bor <- object$model[nrow(object$model)/11*0.8, dep_var]
    # dataset for prediction
    values <- tibble(
        # hier muss die ref class noch einmal anders eingelesen werden 
        # andere möglichkeit nrow()/11*0.8
        #Ref =  new_data$p_value_E_G,
        PRED = object$fitted.values,
        dependent = object$model[, dep_var])
    
    # max and min values of the response variable to correct for it 
    dep_max <- max(values$dependent)
    dep_min <- min(values$dependent)
    
    # computing corrected predictions
    values <- values%>%
        dplyr::mutate(PRED_cor = case_when(
            PRED <= dep_min ~ dep_min + 1e-12,
            PRED >= dep_max ~ dep_max - 1e-12,
            TRUE ~ PRED))
    
    # RMSE and corrected RMSE on the full dataset
    RMSE <- sqrt((sum((values$PRED - values$dependent)^2)/nrow(values)))
    RMSE_cor <- sqrt((sum((values$PRED_cor - values$dependent)^2)/nrow(values)))
    
    # smaller dataset only the intersting part, correction is at 
    values_0.2 <- values %>%
        dplyr::filter(dependent >= dep_var_int_bor)
    RMSE_0.2 <- sqrt((sum((values_0.2$PRED - values_0.2$dependent)^2)/nrow(values_0.2)))
    RMSE_cor_0.2 <- sqrt((sum((values_0.2$PRED_cor - values_0.2$dependent)^2)/nrow(values_0.2)))
    
    # save call and case
    mod_call <- object$call$form[3]
    case <- object$call$data
    
    # values for the plots 
    values <- values%>%
        dplyr::filter(dependent %in% seq(0.001, 1.0 , 0.001))
    
    mod_sum <- tibble(model = as.character(mod_call),
                      RMSE = as.numeric(RMSE),
                      RMSE_cor = as.numeric(RMSE_cor),
                      RMSE_0.2 = as.numeric(RMSE_0.2), 
                      RMSE_cor_0.2 = as.numeric(RMSE_cor_0.2), 
                      case = as.character(case),
                      obs = list(values %>% dplyr::select(dependent)),
                      pred = list(values %>% dplyr::select(PRED)),
                      pred_cor = list(values %>% dplyr::select(PRED_cor))
    )
    return(mod_sum)
}

bind_model_metrics <- function(new_metrics, old_metrics =  model_metrics_E_J) {
    model_metrics_E_J <- rbind(old_metrics,
                               new_metrics)
    
    return(model_metrics_E_J)
}

# Delete Function 
delete_fun <- function(){
    rm(list=setdiff(ls(.GlobalEnv), important_things), envir = .GlobalEnv)
}


#-----------------------
# Preliminary 
#-----------------------
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

# - Boxcox Transformation ----

# case_1 
lambda_x_case_1_stat <- Rfast::bc(data_case_1$stat_Fisher_E_J)
lambda_x_case_1_p <- Rfast::bc(data_case_1$p_value_Fisher_E_J)

data_case_1 <- data_case_1 %>% mutate(
    stat_Fisher_E_J_bc = ((stat_Fisher_E_J^lambda_x_case_1_stat)-1)/lambda_x_case_1_stat,
    p_value_Fisher_E_J_bc = ((p_value_Fisher_E_J^lambda_x_case_1_p)-1)/lambda_x_case_1_p
    )
    

# case_2 
lambda_x_case_2_stat <- Rfast::bc(data_case_2$stat_Fisher_E_J)
lambda_x_case_2_p <- Rfast::bc(data_case_2$p_value_Fisher_E_J)

data_case_2 <- data_case_2 %>% mutate(
    stat_Fisher_E_J_bc = ((stat_Fisher_E_J^lambda_x_case_2_stat)-1)/lambda_x_case_2_stat,
    p_value_Fisher_E_J_bc = ((p_value_Fisher_E_J^lambda_x_case_2_p)-1)/lambda_x_case_2_p
)

# case_3 
lambda_x_case_3_stat <- Rfast::bc(data_case_3$stat_Fisher_E_J)
lambda_x_case_3_p <- Rfast::bc(data_case_3$p_value_Fisher_E_J)

data_case_3 <- data_case_3 %>% mutate(
    stat_Fisher_E_J_bc = ((stat_Fisher_E_J^lambda_x_case_3_stat)-1)/lambda_x_case_3_stat,
    p_value_Fisher_E_J_bc = ((p_value_Fisher_E_J^lambda_x_case_3_p)-1)/lambda_x_case_3_p
)


#-----------------------
# Set up metrics
#-----------------------
model_metrics_E_J <- NULL

# save important things
important_things <- c(ls(), 
                      "important_things")


#-----------------------
# Models - Case 1
#-----------------------
tictoc::tic()



models_poly <- list(
    # functional form: poly(t, p) + (1/k)
    mod_E_J_case.1_p_3 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 3) + I(1/k),
                             data = data_case_1),

    mod_E_J_case.1_p_4 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 4) + I(1/k),
                             data = data_case_1),

    mod_E_J_case.1_p_5 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 5) + I(1/k),
                             data = data_case_1),

    mod_E_J_case.1_p_6 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 6) + I(1/k),
                             data = data_case_1), 

    # functional form: poly(t, p) + (1/k) + poly(t, p)*(1/k)

    mod_E_J_case.1_p_3_3 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 3) + I(1/k) + poly(k, 3)*(1/k),
                               data = data_case_1),

    mod_E_J_case.1_p_4_4 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 4) + I(1/k) + poly(k, 4)*(1/k),
                               data = data_case_1),

    mod_E_J_case.1_p_5_5 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 5) + I(1/k) + poly(k, 5)*(1/k),
                               data = data_case_1),

    mod_E_J_case.1_p_6_6 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 6) + I(1/k) + poly(k, 6)*(1/k),
                               data = data_case_1)
)

# save metrics
model_metrics_E_J <- bind_model_metrics(models_poly %>% purrr::map_dfr(metric_fun))


# delete model
delete_fun()

models_poly_2 <- list(
    # functional form: poly(t, p) + log(k) + poly(t, p)*log(k)
    mod_E_J_case.1_p_3_log <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 3) + log(k) + poly(k, 3)*log(k),
                                 data = data_case_1),

    mod_E_J_case.1_p_4_log <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 4) + log(k) + poly(k, 4)*log(k),
                                 data = data_case_1),

    mod_E_J_case.1_p_5_log <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 5) + log(k) + poly(k, 5)*log(k),
                                 data = data_case_1),

    mod_E_J_case.1_p_6_log <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 6) + log(k) + poly(k, 6)*log(k),
                                 data = data_case_1),

    # functional form: poly(t, p) + k + (1/k)
    mod_E_J_case.1_p_3_k_1 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 3) + k + I(1/k),
                                 data = data_case_1),

    mod_E_J_case.1_p_4_k_1 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 4) + k + I(1/k),
                                 data = data_case_1),

    mod_E_J_case.1_p_5_k_1 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 5) + k + I(1/k),
                                 data = data_case_1),

    mod_E_J_case.1_p_6_k_1 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 6) + k + I(1/k),
                                 data = data_case_1)
)

# save metrics
model_metrics_E_J <- bind_model_metrics(models_poly_2 %>% purrr::map_dfr(metric_fun))

# delete model
delete_fun()

# log()
models_log <- list(
    mod_E_J_case.1_log_1 = lm(p_value_Fisher_E_J ~ log(stat_Fisher_E_J) + k,
                             data = data_case_1),
    
    mod_E_J_case.1_log_2 = lm(p_value_Fisher_E_J ~ log(stat_Fisher_E_J) * k,
                           data = data_case_1),
    
    mod_E_J_case.1_poly_log_3 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 3) + log(k),
                           data = data_case_1),
    
    mod_E_J_case.1_poly_log_4 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 4) + log(k),
                                data = data_case_1),
    
    mod_E_J_case.1_poly_log_5 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 5) + log(k),
                                data = data_case_1),
    
    mod_E_J_case.1_poly_log_6 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 6) + log(k),
                                data = data_case_1),
    
    mod_E_J_case.1_poly_log_7 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 7) + log(k),
                                data = data_case_1),
    
    mod_E_J_case.1_poly_log_8 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 8) + log(k),
                                      data = data_case_1), 
    
    mod_E_J_case.1_poly_log_9 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 9) + log(k),
                                      data = data_case_1),
    
    mod_E_J_case.1_poly_log_10 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 10) + log(k),
                                       data = data_case_1),
    
    mod_E_J_case.1_poly_log_m_3 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 3) * log(k),
                                data = data_case_1),
    
    mod_E_J_case.1_poly_log_m_4 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 4) * log(k),
                                data = data_case_1),
    
    mod_E_J_case.1_poly_log_m_5 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 5) * log(k),
                                data = data_case_1),
    
    mod_E_J_case.1_poly_log_m_6 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 6)  *log(k),
                                data = data_case_1),
    
    mod_E_J_case.1_poly_log_m_7 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 7) * log(k),
                                data = data_case_1),
    
    mod_E_J_case.1_poly_log_m_8 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 8) * log(k),
                                      data = data_case_1), 
    
    mod_E_J_case.1_poly_log_m_9 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 9) * log(k),
                                      data = data_case_1),
    
    mod_E_J_case.1_poly_log_m_10 <- lm(p_value_Fisher_E_J ~ poly(log(stat_Fisher_E_J), 10) * log(k),
                                      data = data_case_1)
)

# save metrics
model_metrics_E_J <- bind_model_metrics(models_log %>% purrr::map_dfr(metric_fun))

# delete model
delete_fun()


models_bc <- list(
    mod_E_J_case.1_bc_poly_log_m_10 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J_bc, 10) * log(k),
                                         data = data_case_1),
    mod_E_J_case.1_bc_poly_m_10 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J_bc, 10) * k, 
                                     data = data_case_1),
    mod_E_J_case.1_bc_poly_10 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J_bc, 10) + k, 
                                     data = data_case_1),
    mod_E_J_case.1_bc_poly_log_10 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J_bc, 10) + log(k), 
                                   data = data_case_1),
    mod_E_J_case.1_bc_poly_log_m_10_I <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J_bc, 10) * log(k) + I(1/k),
                                         data = data_case_1),
    mod_E_J_case.1_bc_poly_m_10_I <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J_bc, 10) * k + I(1/k),
                                           data = data_case_1),
    mod_E_J_case.1_bc_poly_log_m_10_poly_m_I <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J_bc, 10) * log(k) + poly(stat_Fisher_E_J_bc, 10)*I(1/k),
                                           data = data_case_1),
    mod_E_J_case.1_bc_poly_log_m_10_poly_m_I_sqrt <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J_bc, 10) * log(k) + poly(stat_Fisher_E_J_bc, 10)*I(1/k) + sqrt(k),
                                                  data = data_case_1),
    mod_E_J_case.1_bc_poly_log_m_10_poly_m_sqrt <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J_bc, 10) * log(k) + poly(stat_Fisher_E_J_bc, 10)*sqrt(k),
                                                  data = data_case_1)
)

model_metrics_E_J <- bind_model_metrics(models_bc %>% purrr::map_dfr(metric_fun))


#--- GAM ---

# Generalized Additive Models
models_gam <- list(
    
    # functional form: poly(t, p) + (1/k)
    gam_case.1_ns_3 <-  mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 3) + ns(I(1/k), 3),
                                     data = data_case_1),
    
    gam_case.1_ns_5 <-  mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 5) + ns(I(1/k), 5),
                                  data = data_case_1),
    
    gam_case.1_ns_7 <-  mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 7) + ns(I(1/k), 7),
                                  data = data_case_1),
    
    gam_case.1_ns_9 <-  mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 9) + ns(I(1/k), 9),
                                  data = data_case_1),
    
    # functional form: poly(t, p) + (1/k) + poly(t, p)*(1/k)
    gam_case.1_ns_3_3 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 3) + I(1/k) + ns(k, 3)*(1/k),
                               data = data_case_1),
    
    gam_case.1_ns_5_5 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 5) + I(1/k) + ns(k, 5)*(1/k),
                                   data = data_case_1),
    
    gam_case.1_ns_7_7 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 7) + I(1/k) + ns(k, 7)*(1/k),
                                   data = data_case_1),
    
    gam_case.1_ns_9_9 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 9) + I(1/k) + ns(k, 9)*(1/k),
                                   data = data_case_1),
    
    # functional form: poly(t, p) + log(k) + poly(t, p)*log(k)
    gam_case.1_ns_3_log <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 3) + log(k) + ns(k, 3)*log(k),
                                     data = data_case_1),
    
    gam_case.1_ns_5_log <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 5) + log(k) + ns(k, 5)*log(k),
                                     data = data_case_1),
    
    gam_case.1_ns_7_log <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 7) + log(k) + ns(k, 7)*log(k),
                                     data = data_case_1),
    
    gam_case.1_ns_9_log <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 9) + log(k) + ns(k, 9)*log(k),
                                     data = data_case_1),
    
    # functional form: poly(t, p) + k + (1/k)
    gam_case.1_p_3_k_1 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 3) + k + I(1/k),
                                 data = data_case_1),
    
    gam_case.1_p_5_k_1 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 5) + k + I(1/k),
                                    data = data_case_1),
    
    gam_case.1_p_7_k_1 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 7) + k + I(1/k),
                                    data = data_case_1),
    
    gam_case.1_p_9_k_1 <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J, 9) + k + I(1/k),
                                    data = data_case_1),
    
    # functional form: ns(log(t, p)) + log(k) 
    gam_case.1_ns_3_log_ad <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 3) + log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_5_log_ad <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 5) + log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_7_log_ad <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 7) + log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_9_log_ad <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 9) + log(k),
                                        data = data_case_1),
    
    # functional form: ns(log(t, p)) * log(k) + (1/k)
    gam_case.1_ns_3_log_mu <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 3) * log(k),
                                     data = data_case_1),
    
    gam_case.1_ns_5_log_mu <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 5) * log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_7_log_mu <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 7) * log(k),
                                        data = data_case_1),
    
    gam_case.1_ns_9_log_mu <- mgcv::gam(p_value_Fisher_E_J ~ ns(log(stat_Fisher_E_J), 9) * log(k),
                                        data = data_case_1),
    
    # bc 
    gam_case.1_bc_poly_log_m_10_poly_m_sqrt <- mgcv::gam(p_value_Fisher_E_J ~ ns(stat_Fisher_E_J_bc, 9) * log(k) + 
                                                             ns(stat_Fisher_E_J_bc, 9)*sqrt(k), data = data_case_1),
    
    # by ~ bc 
    gam_case.1_bc_poly_log_m_10_poly_m_sqrt <- mgcv::gam(p_value_Fisher_E_J_bc ~ ns(stat_Fisher_E_J_bc, 3) * log(k) + 
                                                                 ns(stat_Fisher_E_J_bc, 3)*sqrt(k), data = data_case_1)

    
)


# save metrics
model_metrics_E_J <- bind_model_metrics(models_gam %>% purrr::map_dfr(metric_fun))


delete_fun()
tictoc::toc()


#-----------------------
# Models - Case 2
#-----------------------

#-----------------------
# Models - Case 3
#-----------------------



#---- saving metrics ----


save(model_metrics_E_J, file = here::here("09_simulation_and_approximation-cdt/model_metrics_E_J.Rdata"))

