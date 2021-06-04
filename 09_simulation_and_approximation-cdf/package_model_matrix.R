#---- Preliminary ---- 
# packages  
source(here::here('01_code/packages/packages.R'))

# final models and lambda
load(here::here('09_simulation_and_approximation-cdf/lambda_package.RData'))
load(here::here('09_simulation_and_approximation-cdf/final_models.RData'))

## wide format
lambda_values %<>%
    tidyr::pivot_wider(names_from = c('side')) %>%
    dplyr::rename(test.type = test) %>%
    dplyr::mutate(test.type = dplyr::recode(test.type, 'e_j' = 'E_J'))

models <- final_models %>%
    dplyr::left_join(lambda_values, by = c('case', 'test.type'))

save(models, file = here::here('09_simulation_and_approximation-cdf/models.RData'))

