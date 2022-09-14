#---- Preliminary ---- 
# packages
source(here::here('01_code/packages/packages.R'))

# load functions 
source(here::here('09_simulation_and_approximation-cdf/estimation_functions.R'))

# parallel
plan(multisession, workers = 3)
# set off warning for random numbers
options(future.rng.onMisuse = "ignore", 
        future.globals.maxSize = 2.147e+9)

model_cleaner = function(object) {
    object$y = c()
    object$model = c()
    
    object$residuals = c()
    object$fitted.values = c()
    object$effects = c()
    object$qr$qr = c()
    object$linear.predictors = c()
    object$weights = c()
    object$prior.weights = c()
    object$data = c()
    # new 
    object$df.residual = c()
    return(object)
}

#-- tibble with models and data ----

# Blabla, ich lad die Scheiße jetzt manuell :D
load("~/Desktop/data_cases.RData")

load("~/Documents/GitHub/cointegration_intern/09_simulation_and_approximation-cdf/lambda_package.RData")

#---- Models ----
library(magrittr)
# E_J
model_EJ_1 <- lm(p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, 13) * k_dummy, 
                 data = data_case_1)
model_EJ_1 %<>% model_cleaner()

model_EJ_2 <- lm(p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, 13) * k_dummy, 
                 data = data_case_2)
model_EJ_2 %<>% model_cleaner()

model_EJ_3 <- lm(p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, 13) * k_dummy, 
                 data = data_case_3)
model_EJ_3 %<>% model_cleaner()

# all
model_all_1 <- lm(p_value_Fisher ~ poly(stat_Fisher_all_bc, 13) * k_dummy, 
                  data = data_case_1)
model_all_1 %<>% model_cleaner()

model_all_2 <- lm(p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, 12) * k_dummy, 
                  data = data_case_2)
model_all_2 %<>% model_cleaner()

model_all_3 <- lm(p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, 13) * k_dummy, 
                  data = data_case_3)
model_all_3 %<>% model_cleaner()


#save(model_EJ_1, model_EJ_2, 
#     file = here::here('09_simulation_and_approximation-cdf/models_lm.RData'))

#load(file = here::here('09_simulation_and_approximation-cdf/models_lm_small.RData'))
save(model_EJ_1, model_EJ_2, model_EJ_3,
     model_all_1, model_all_2, model_all_3,
     file = here::here('09_simulation_and_approximation-cdf/models_lm_small.RData'))


