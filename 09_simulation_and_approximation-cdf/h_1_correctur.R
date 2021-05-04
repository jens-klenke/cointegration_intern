# packages 
source(here::here('01_code/packages/packages.R'))

# functions 
source(here::here('09_simulation_and_approximation-cdf/func_p_approximation.R'))

# load models 
load(here::here('09_simulation_and_approximation-cdf/models_package.RData'))

# H_1 limit 
load(here::here('09_simulation_and_approximation-cdf/H_1-values.RData'))

# functions 
ciritical_fun <- function(case_s, test_s){
    crit_val %>%
        dplyr::filter(case == case_s, 
                      test == test_s) %>%
        dplyr::ungroup() %>%
        dplyr::select(-c(case, test))
}

# storing all critical_val in one tibble 

crit_val <- dplyr::bind_rows(
    crit_val_all_10 %>%
        dplyr::rename(crit_val = 'stat_Fisher_all') %>%
        dplyr::mutate(test = 'all'),
    
    crit_val_e_j_10 %>%
        dplyr::mutate(test = 'e_j') %>%
        dplyr::rename(crit_val = 'stat_Fisher_E_J'))


# adding critical val to tibble

models %<>%
    dplyr::mutate(critical = purrr::map2(case, test, ciritical_fun))




