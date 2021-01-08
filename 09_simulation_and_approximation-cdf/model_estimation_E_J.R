### functions 
metric_fun <- function(object){
    model <- substitute(object)
    
    values <- tibble(PRED = object$fitted.values,
                        dependent = object$model[,'p_value_Fisher_E_J'])%>%
        dplyr::mutate(PRED_cor = case_when(
            PRED <= 0 ~ 1e-12,
            PRED >= 1 ~ 1-1e-12,
            TRUE ~ PRED))
    
    
    RMSE <- sqrt((sum((values$PRED - values$dependent)^2)/nrow(values)))
    RMSE_cor <- sqrt((sum((values$PRED_cor - values$dependent)^2)/nrow(values)))
    
    values_0.2 <- values %>%
        dplyr::filter(dependent >= 0.8)
    
    RMSE_0.2 <- sqrt((sum((values_0.2$PRED - values_0.2$dependent)^2)/nrow(values_0.2)))
    RMSE_cor_0.2 <- sqrt((sum((values_0.2$PRED_cor - values_0.2$dependent)^2)/nrow(values_0.2)))
    
    mod_call <- object$call$form[3]
    case <- object$call$data
    
    metric <- c(as.character(model), as.numeric(RMSE), as.numeric(RMSE_cor), 
                as.numeric(RMSE_0.2), as.numeric(RMSE_cor_0.2), as.character(mod_call), as.character(case))
    
    return(metric)
}

## packages
source(here::here('01_code/packages/packages.R'))

## load simulation data 

# jens 
if(Sys.info()['nodename'] == "DESKTOP-ATT92OH"){
    load('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_1_m.RData')
} else if(Sys.info()['nodename'] == "OEK-TS01"){ # Server
    load('D:\\Klenke\\Data_1_m.RData')
}

# split dataset in cases
data_case_1 <- Data %>%
    dplyr::filter(case == 1)

data_case_2 <- Data %>%
    dplyr::filter(case == 2)

data_case_3 <- Data %>%
    dplyr::filter(case == 3)

### set up metrics
model_metrics_E_J <- NULL

tictoc::tic()

#### models  case 1 ####
# functional form: poly(t, p) + (1/k)


mod_E_J_case.1_p_3 = lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 3) + I(1/k),
                                   data = data_case_1)


metric_fun(mod_E_J_case.1_p_3)

mod_E_J_case.1_p_3 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 3) + I(1/k),
                         data = data_case_1)

mod_E_J_case.1_p_4 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 4) + I(1/k),
                                   data = data_case_1)

mod_E_J_case.1_p_5 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 5) + I(1/k),
                                   data = data_case_1)
mod_E_J_case.1_p_6 <- lm(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 6) + I(1/k),
                                   data = data_case_1)

# save metrics
model_metrics_E_J <- rbind(model_metrics ,
                       metric_fun(mod_E_J_case.1_p_3, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_4, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_5, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_6, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'))

# delete model
rm(list = c('mod_E_J_case.1_p_3', 'mod_E_J_case.1_p_4', 'mod_E_J_case.1_p_5', 'mod_E_J_case.1_p_6'))


# functional form: poly(t, p) + (1/k) + poly(t, p)*(1/k)
mod_E_J_case.1_p_3_3 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 3) + I(1/k) + poly(k, 3)*(1/k),
                                     data = data_case_1, 
                                     method = 'lm', 
                                     trControl = CV_control)

mod_E_J_case.1_p_4_4 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 4) + I(1/k) + poly(k, 4)*(1/k),
                                     data = data_case_1, 
                                     method = 'lm', 
                                     trControl = CV_control)

mod_E_J_case.1_p_5_5 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 5) + I(1/k) + poly(k, 5)*(1/k),
                                     data = data_case_1, 
                                     method = 'lm', 
                                     trControl = CV_control)

mod_E_J_case.1_p_6_6 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 6) + I(1/k) + poly(k, 6)*(1/k),
                                     data = data_case_1, 
                                     method = 'lm', 
                                     trControl = CV_control)

# save metrics
model_metrics <- rbind(model_metrics ,
                       metric_fun(mod_E_J_case.1_p_3_3, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_4_4, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_5_5, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_6_6, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'))

# delete model
rm(list = c('mod_E_J_case.1_p_3_3', 'mod_E_J_case.1_p_4_4', 'mod_E_J_case.1_p_5_5', 'mod_E_J_case.1_p_6_6'))


####  chi^p ln  
mod_E_J_case.1_p_3_log <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 3) + log(k) + poly(k, 3)*log(k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control)

mod_E_J_case.1_p_4_log <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 4) + log(k) + poly(k, 4)*log(k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control)

mod_E_J_case.1_p_5_log <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 5) + log(k) + poly(k, 5)*log(k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control)

mod_E_J_case.1_p_6_log <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 6) + log(k) + poly(k, 6)*log(k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control)

# save metrics
model_metrics <- rbind(model_metrics ,
                       metric_fun(mod_E_J_case.1_p_3_log, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_4_log, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_5_log, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_6_log, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'))

# delete model
rm(list = c('mod_E_J_case.1_p_3_log', 'mod_E_J_case.1_p_4_log', 'mod_E_J_case.1_p_5_log', 'mod_E_J_case.1_p_6_log'))

####  chi^p + k + (1/k)  
mod_E_J_case.1_p_3_k_1 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 3) + k + I(1/k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control)

mod_E_J_case.1_p_4_k_1 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 4) + k + I(1/k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control)

mod_E_J_case.1_p_5_k_1 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 5) + k + I(1/k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control)

mod_E_J_case.1_p_6_k_1 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 6) + k + I(1/k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control)

# save metrics
model_metrics <- rbind(model_metrics ,
                       metric_fun(mod_E_J_case.1_p_3_k_1, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_4_k_1, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_5_k_1, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'),
                       metric_fun(mod_E_J_case.1_p_6_k_1, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J'))

# delete model
rm(list = c('mod_E_J_case.1_p_3_k_1', 'mod_E_J_case.1_p_4_k_1', 'mod_E_J_case.1_p_5_k_1', 'mod_E_J_case.1_p_6_k_1'))


colnames(model_metrics) <- c('model', 'RMSE', 'RMSE_p<0.2', 'formula')

### GAM 
mod_E_G_case.1_gam_3 <- mgcv::gam(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 6) + ns(k, 6) + ns(stat_E_G*k, 6),
                                  data = data_case_1)
    
    
    
metric_fun(mod_E_G_case.1_gam_3, test_E_J_data_case_1, dep_VAR = 'p_value_Fisher_E_J')


tictoc::toc()




############################ Case 2 ##########################

fit_case_2 <- lm(p_value_E_G ~ poly(stat_E_G, 3) + poly(k, 3) ,  data = data_case_2)
summary(fit_case_2)


############################ Case 3 ##########################

fit_case_3 <- lm(p_value_E_G ~ poly(stat_E_G, 3) + poly(k, 3) ,  data = data_case_3)
summary(fit_case_3)




