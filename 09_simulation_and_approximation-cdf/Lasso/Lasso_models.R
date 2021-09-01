#load("/Users/Janine/Desktop/BayerHanck/data_cases.RData")

library(glmnet)
#---- EJ models ----

# case 1
x_EJ_1 <- model.matrix(p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, 13) * k_dummy, 
                       data = data_case_1)[, -1]
y_bc <- data_case_1$p_value_Fisher_bc

lasso_EJ_1 <- glmnet::glmnet(x = x_EJ_1, y = y_bc, alpha = 1)

# case 2
x_EJ_2 <- model.matrix(p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, 13) * k_dummy, 
                       data = data_case_2)[, -1]

lasso_EJ_2 <- glmnet::glmnet(x = x_EJ_2, y = y_bc, alpha = 1)

# case 3
x_EJ_3 <- model.matrix(p_value_Fisher_bc ~ poly(stat_Fisher_E_J_bc, 13) * k_dummy, 
                       data = data_case_3)[, -1]

lasso_EJ_3 <- glmnet::glmnet(x = x_EJ_3, y = y_bc, alpha = 1)


#---- all models ----

# case 1
x_all_1 <- model.matrix(p_value_Fisher ~ poly(stat_Fisher_all_bc, 13) * k_dummy, 
                       data = data_case_1)[, -1]
y <- data_case_1$p_value_Fisher

lasso_all_1 <- glmnet::glmnet(x = x_all_1, y = y, alpha = 1)

# case 2
x_all_2 <- model.matrix(p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, 12) * k_dummy, 
                       data = data_case_2)[, -1]

lasso_all_2 <- glmnet::glmnet(x = x_all_2, y = y_bc, alpha = 1)

# case 3
x_all_3 <- model.matrix(p_value_Fisher_bc ~ poly(stat_Fisher_all_bc, 13) * k_dummy, 
                       data = data_case_3)[, -1]

lasso_all_3 <- glmnet::glmnet(x = x_all_3, y = y_bc, alpha = 1)


# ---- save ----

save(lasso_EJ_1, lasso_EJ_2, lasso_EJ_3, lasso_all_1, lasso_all_2, lasso_all_3,
     file = here::here("09_simulation_and_approximation-cdf/Lasso/models.RData"))

save(y_bc, y, x_EJ_1, x_EJ_2, x_EJ_3, x_all_1, x_all_2, x_all_3,
     file = here::here("09_simulation_and_approximation-cdf/Lasso/model_matrix.RData"))
