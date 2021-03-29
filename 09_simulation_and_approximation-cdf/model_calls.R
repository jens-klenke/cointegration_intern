#---- Preliminary ---- 
source(here::here('01_code/packages/packages.R'))

# load metrics
load(here::here("09_simulation_and_approximation-cdf/model_metrics_ALL_server.Rdata"))

load(here::here("09_simulation_and_approximation-cdf/model_metrics_E_J_server.Rdata"))

# table functional form 
model_functions <- model_metrics_ALL %>%
    dplyr::filter(case == 'data_case_1') %>%
    dplyr::select(model)

# final models 

A <- model_metrics_ALL %>%
    dplyr::filter(case == 'data_case_1') %>%
    dplyr::slice_min(RMSE_cor_0.2, n = 5) %>%
    dplyr::select(-c(pred, case)) 

A_2 <- model_metrics_ALL %>%
    dplyr::filter(case == 'data_case_2') %>%
    dplyr::slice_min(RMSE_cor_0.2, n = 5) %>%
    dplyr::select(-c(pred, case))

A_3 <- model_metrics_ALL %>%
    dplyr::filter(case == 'data_case_3') %>%
    dplyr::slice_min(RMSE_cor_0.2, n = 5) %>%
    dplyr::select(-c(pred, case))

EG_J <- model_metrics_E_J %>%
    dplyr::filter(case == 'data_case_1') %>%
    dplyr::slice_min(RMSE_cor_0.2, n = 5) %>%
    dplyr::select(-c(pred, case))





sum_var <- papeR::summarise( as.data.frame(fifa_data%>%dplyr::select(value_eur, wage_eur, overall, year, age, potential )), 
                             group = 'year', p_value = FALSE, quantil = FALSE )

sum_var <- sum_var[,1:((dim(sum_var)[2])-2)]

sum_var%>%kable()%>%kable_styling()

#knitr::kable(sum_var, digits = 2, caption = '\\label{tab:sum} Summary of some important variables for the 2019 FIFA edition', format.args = list(decimal.mark = '.', big.mark = " "), col.names = c('', 'year', '', 'N', ' ', 'mean', 'sd' ), row.names = FALSE,  booktabs = TRUE, linesep = "")%>%
# kable_styling(latex_options = c("striped", "hold_position"))

