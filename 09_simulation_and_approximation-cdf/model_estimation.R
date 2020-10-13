## packages
source(here::here('01_code/packages/packages.R'))

## load simulation data 
Data <- base::readRDS(here::here('09_simulation_and_approximation-cdf/Data.rds'))

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

tictoc::tic()
#### models  case 1 ####
# functional form: poly(t, p) + (1/k)
mod_E_G_case.1_p_3 <- caret::train(p_value_E_G ~ poly(stat_E_G, 3) + I(1/k),
                               data = data_case_1, 
                               method = 'lm', 
                               trControl = CV_control)

mod_E_G_case.1_p_3$results
mod_E_G_case.1_p_3$call$form

mod_E_G_case.1_p_4 <- caret::train(p_value_E_G ~ poly(stat_E_G, 4) + I(1/k),
                               data = data_case_1, 
                               method = 'lm', 
                               trControl = CV_control )

mod_E_G_case.1_p_5 <- caret::train(p_value_E_G ~ poly(stat_E_G, 5) + I(1/k),
                               data = data_case_1, 
                               method = 'lm', 
                               trControl = CV_control )

mod_E_G_case.1_p_6 <- caret::train(p_value_E_G ~ poly(stat_E_G, 6) + I(1/k),
                               data = data_case_1, 
                               method = 'lm', 
                               trControl = CV_control )

# functional form: poly(t, p) + (1/k) + poly(t, p)*(1/k)
mod_E_G_case.1_p_3_3 <- caret::train(p_value_E_G ~ poly(stat_E_G, 3) + I(1/k) + poly(k, 3)*(1/k),
                               data = data_case_1, 
                               method = 'lm', 
                               trControl = CV_control )

mod_E_G_case.1_p_4_4 <- caret::train(p_value_E_G ~ poly(stat_E_G, 4) + I(1/k) + poly(k, 4)*(1/k),
                               data = data_case_1, 
                               method = 'lm', 
                               trControl = CV_control )

mod_E_G_case.1_p_5_5 <- caret::train(p_value_E_G ~ poly(stat_E_G, 5) + I(1/k) + poly(k, 5)*(1/k),
                               data = data_case_1, 
                               method = 'lm', 
                               trControl = CV_control )

mod_E_G_case.1_p_6_6 <- caret::train(p_value_E_G ~ poly(stat_E_G, 6) + I(1/k) + poly(k, 6)*(1/k),
                               data = data_case_1, 
                               method = 'lm', 
                               trControl = CV_control )


####  chi^p ln  
mod_E_G_case.1_p_3_log <- caret::train(p_value_E_G ~ poly(stat_E_G, 3) + log(k) + poly(k, 3)*log(k),
                                 data = data_case_1, 
                                 method = 'lm', 
                                 trControl = CV_control )

mod_E_G_case.1_p_4_log <- caret::train(p_value_E_G ~ poly(stat_E_G, 4) + log(k) + poly(k, 4)*log(k),
                                 data = data_case_1, 
                                 method = 'lm', 
                                 trControl = CV_control )

mod_E_G_case.1_p_5_log <- caret::train(p_value_E_G ~ poly(stat_E_G, 5) + log(k) + poly(k, 5)*log(k),
                                 data = data_case_1, 
                                 method = 'lm', 
                                 trControl = CV_control )

mod_E_G_case.1_p_6_log <- caret::train(p_value_E_G ~ poly(stat_E_G, 6) + log(k) + poly(k, 6)*log(k),
                                 data = data_case_1, 
                                 method = 'lm', 
                                 trControl = CV_control )

####  chi^p + k + (1/k)  
mod_E_G_case.1_p_3_k_1 <- caret::train(p_value_E_G ~ poly(stat_E_G, 3) + k + I(1/k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control )

mod_E_G_case.1_p_4_k_1 <- caret::train(p_value_E_G ~ poly(stat_E_G, 4) + k + I(1/k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control )

mod_E_G_case.1_p_5_k_1 <- caret::train(p_value_E_G ~ poly(stat_E_G, 5) + k + I(1/k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control )

mod_E_G_case.1_p_6_k_1 <- caret::train(p_value_E_G ~ poly(stat_E_G, 6) + k + I(1/k),
                                       data = data_case_1, 
                                       method = 'lm', 
                                       trControl = CV_control )

tictoc::toc()


## case 2

fit_case_2 <- lm(p_value_E_G ~ poly(stat_E_G, 3) + poly(k, 3) ,  data = data_case_2)
summary(fit_case_2)

fit_case_3 <- lm(p_value_E_G ~ poly(stat_E_G, 3) + poly(k, 3) ,  data = data_case_3)
summary(fit_case_3)



