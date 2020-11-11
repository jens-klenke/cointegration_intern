## packages
source(here::here('01_code/packages/packages.R'))


load(here::here('09_simulation_and_approximation-cdf/sysdata.rda'))


own_sim <- readRDS(here::here('09_simulation_and_approximation-cdf/Data.Rds'))


own_sim_critical_val <- own_sim %>% 
    dplyr::filter(p_value_Fisher_all == 0.95)%>%
    dplyr::select(p_value_Fisher_all, stat_Fisher_all,  k, case)%>%
    dplyr::arrange(case, k)%>%
    dplyr::mutate(c_val = c(crit_val_2_0.05))%>%
    dplyr::mutate(diff = stat_Fisher_all - c_val, 
                  diff_procent = diff/c_val * 100)




