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

models$model[[3]]$formula


#save(model_EJ_1, model_EJ_2, 
#     file = here::here('09_simulation_and_approximation-cdf/models_lm.RData'))

save(model_EJ_1, model_EJ_2, 
     file = here::here('09_simulation_and_approximation-cdf/models_lm_small.RData'))
