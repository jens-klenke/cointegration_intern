load(here::here("09_simulation_and_approximation-cdf/Lasso/lasso_models.RData"))
load("~/Desktop/Lasso/model_matrix.RData")
load(here::here("09_simulation_and_approximation-cdf/Lasso/lambda_p.RData"))

library(glmnet)
library(magrittr)
library(dplyr)

# Functions
lasso_eval <- function(mod, x_matrix, y = y) {
  dep_var <- mod$call$y
  lambda <- mod$lambda %>% 
    dplyr::last()
  fitted_values <- glmnet::predict.glmnet(mod, newx = x_matrix, s = lambda)
  
  
  values <- tibble(
    PRED = if(stringr::str_detect(dep_var, '_bc')){ 
      invBoxCox(fitted_values)
    } else if(stringr::str_detect(dep_var, '_lg')){
      exp(fitted_values)
    } else {fitted_values[, 1]}, # y_hat 
    dependent = y
  )
  
  # computing corrected predictions
  values %<>% dplyr::mutate(PRED_cor = case_when(
    PRED <= 0 ~ 1e-12,
    PRED >= 1 ~ 1 - 1e-12,
    TRUE ~ PRED))
  
  # RMSE and corrected RMSE on full dataset
  RMSE <- sqrt(sum((values$PRED - values$dependent)^2)/nrow(values))
  RMSE_cor <- sqrt(sum((values$PRED_cor - values$dependent)^2)/nrow(values))
  
  # subset dataset: only the interesting part
  values_0.2 <- values %>%
    dplyr::filter(dependent >= 0.2)
  RMSE_0.2 <- sqrt(sum((values_0.2$PRED - values_0.2$dependent)^2)/nrow(values_0.2))
  RMSE_cor_0.2 <- sqrt(sum((values_0.2$PRED_cor - values_0.2$dependent)^2)/nrow(values_0.2))
  
  tibble(model  = list(mod),
         RMSE = RMSE,
         RMSE_cor = RMSE_cor,
         RMSE_0.2 = RMSE_0.2, 
         RMSE_cor_0.2 = RMSE_cor_0.2
  )
}
invBoxCox <- function(x){
  x <- if (lambda_p == 0) exp(as.complex(x)) else (lambda_p*as.complex(x) + 1)^(1/lambda_p)
  return(Re(x))
}

#---- all ----
pred_all_1 <- lasso_eval(mod = lasso_all_1, x_matrix = x_all_1, y = y)
pred_all_2 <- lasso_eval(mod = lasso_all_2, x_matrix = x_all_2, y = y)
pred_all_3 <- lasso_eval(mod = lasso_all_3, x_matrix = x_all_3, y = y)

#---- EJ ----
pred_EJ_1 <- lasso_eval(mod = lasso_EJ_1, x_matrix = x_EJ_1, y = y)
pred_EJ_2 <- lasso_eval(mod = lasso_EJ_2, x_matrix = x_EJ_2, y = y)
pred_EJ_3 <- lasso_eval(mod = lasso_EJ_3, x_matrix = x_EJ_3, y = y)


calls <- c(
  "$p = c + \\poly\\left( \\bc(t), 13 \\right) * k\\_d$", 
  "$\\bc(p) = c + \\poly\\left( \\bc(t), 12 \\right) * k\\_d$",
  "$\\bc(p) = c + \\poly\\left( \\bc(t), 13 \\right) * k\\_d$",
  "$\\bc(p) = c + \\poly\\left( \\bc(t), 13 \\right) * k\\_d$",
  "$\\bc(p) = c + \\poly\\left( \\bc(t), 13 \\right) * k\\_d$",
  "$\\bc(p) = c + \\poly\\left( \\bc(t), 13 \\right) * k\\_d$"
)

lasso <- dplyr::bind_rows(pred_all_1, pred_all_2, pred_all_3,
                          pred_EJ_1, pred_EJ_2, pred_EJ_3) %>%
  dplyr::select(-model) %>%
  dplyr::bind_cols(test.type = rep(c("all", "E_J"), each = 3), 
                   case = rep(1:3, 2), calls = calls, .)

save(lasso, file = here::here("09_simulation_and_approximation-cdf/Lasso/lasso_table.RData"))

