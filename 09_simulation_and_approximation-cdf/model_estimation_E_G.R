### functions 
metric_fun <- function(object, new_data_0.2, dep_VAR){
    model <- substitute(object)
    RMSE <- object$results$RMSE
    
    PRED <- predict(object, new_data_0.2)
    RMSE_0.2 <- sqrt((sum((new_data_0.2[, dep_VAR] - PRED)^2)/length(PRED)))
    
    call <- object$call$form[3]
    
    metric <- c(as.character(model), as.numeric(RMSE), as.numeric(RMSE_0.2), as.character(call))
    
    return(metric)
}

## packages
source(here::here('01_code/packages/packages.R'))

## load simulation data 
#Data <- base::readRDS(here::here('09_simulation_and_approximation-cdf/Data.rds'))
if(Sys.info()['user'] == "Jens-"){
    Data <- base::readRDS('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data.rds')
}


# split dataset in cases
data_case_1 <- Data %>%
    dplyr::filter(case == 1)

data_case_2 <- Data %>%
    dplyr::filter(case == 2)

data_case_3 <- Data %>%
    dplyr::filter(case == 3)

### train caret 
CV_control <- trainControl(method = "cv", number = 10) 

#### E_G_model #### 

# test dataset for 0.2

test_E_J_data_case_1 <- Data %>%
    dplyr::filter(case == 1, 
                  p_value_Fisher_E_J <= 0.2)

test_E_J_data_case_2 <- Data %>%
    dplyr::filter(case == 2, 
                  p_value_Fisher_E_J <= 0.2)

test_E_J_data_case_3 <- Data %>%
    dplyr::filter(case == 3, 
                  p_value_Fisher_E_J <= 0.2)

### set up metrics
model_metrics <- NULL

tictoc::tic()

#### models  case 1 ####
# functional form: poly(t, p) + (1/k)
mod_E_J_case.1_p_3 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 3) + I(1/k),
                                   data = data_case_1, 
                                   method = 'lm', 
                                   trControl = CV_control)

mod_E_J_case.1_p_4 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 4) + I(1/k),
                                   data = data_case_1, 
                                   method = 'lm', 
                                   trControl = CV_control)

mod_E_J_case.1_p_5 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 5) + I(1/k),
                                   data = data_case_1, 
                                   method = 'lm', 
                                   trControl = CV_control)

mod_E_J_case.1_p_6 <- caret::train(p_value_Fisher_E_J ~ poly(stat_Fisher_E_J, 6) + I(1/k),
                                   data = data_case_1, 
                                   method = 'lm', 
                                   trControl = CV_control)

# save metrics
model_metrics <- rbind(model_metrics ,
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




