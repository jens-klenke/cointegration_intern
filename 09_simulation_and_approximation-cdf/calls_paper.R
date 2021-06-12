calls <- c('$p = \\poly\\left( \\bc(t), power \\right)$', '$p = \\poly\\left( \\bc(t), power \\right) + k$', 
  '$p = \\poly\\left( \\bc(t), power \\right) * k$', '$p = \\poly\\left( \\bc(t), power \\right) + \\log(k)$',
  '$p = \\poly\\left( \\bc(t), power \\right) * \\log(k)$', '$p = \\poly\\left( \\bc(t), power \\right) + k\\_d$',
  '$p = \\poly\\left( \\bc(t), power \\right) * k\\_d$', '$\\log(p) = \\poly\\left( \\bc(t), power \\right)$',
  '$\\log(p) = \\poly\\left( \\bc(t), power \\right) + k$', '$\\log(p) = \\poly\\left( \\bc(t), power \\right) * k$', 
  '$\\log(p) = \\poly\\left( \\bc(t), power \\right) + \\log(k)$', '$\\log(p) = \\poly\\left( \\bc(t), power \\right) * \\log(k)$',
  '$\\log(p) = \\poly\\left( \\bc(t), power \\right) + k\\_d$', '$\\log(p) = \\poly\\left( \\bc(t), power \\right) * k\\_d$', 
  '$\\bc(p) = \\poly\\left( \\bc(t), power \\right)$', '$\\bc(p) = \\poly\\left( \\bc(t), power \\right) + k$', 
  '$\\bc(p) = \\poly\\left( \\bc(t), power \\right) * k$', '$\\bc(p) = \\poly\\left( \\bc(t), power \\right) + \\log(k)$', 
  '$\\bc(p) = \\poly\\left( \\bc(t), power \\right) * \\log(k)$', '$\\bc(p) = \\poly\\left( \\bc(t), power \\right) + k\\_d$', 
  '$\\bc(p) = \\poly\\left( \\bc(t), power \\right) * k\\_d$')

expo <- 3:13

A <- expand_grid(calls, expo) %>%
    # functional call, merge of power and call
    dplyr::mutate(dplyr::across(calls, str_replace_all, "power", 
                                as.character(.$expo))) %>%
    dplyr::select(calls)

# load metrics
load(here::here('09_simulation_and_approximation-cdf/server_results.RData'))

A %<>%
    dplyr::pull(calls)

table_all_case_1 %<>%
    dplyr::mutate(calls = A)

table_all_case_2 %<>%
    dplyr::mutate(calls = A)

table_all_case_3 %<>%
    dplyr::mutate(calls = A)

table_E_J_case_1 %<>%
    dplyr::mutate(calls = A)

table_E_J_case_2 %<>%
    dplyr::mutate(calls = A)

table_E_J_case_3 %<>%
    dplyr::mutate(calls = A)

save(table_all_case_1, table_all_case_2, table_all_case_3, 
     table_E_J_case_1, table_E_J_case_2, table_E_J_case_3,
    file = here::here('09_simulation_and_approximation-cdf/server_results.RData'))


##

best_5_table_paper <- function(data){
    data%>%
        dplyr::slice_min(RMSE_cor_0.2, n = 5)%>%
        dplyr::select(-c(formula, expo, model))
}

best_5_table_paper(table_all_case_1)


