#------ functions

# inverse BoxCox function
invBoxCox <- function(x){
    x <- if (lambda_p == 0) exp(as.complex(x)) else (lambda_p*as.complex(x) + 1)^(1/lambda_p)
    return(Re(x))
}

# BoxCox + log Transformation
bc_log_fun <- function(data) {
    case <- deparse(substitute(data)) %>% 
        stringr::str_extract("[:digit:]")
    lambda_stat_E_J <- data %>% 
        dplyr::pull(stat_Fisher_E_J) %>%
        Rfast::bc()
    lambda_p <- data %>%
        dplyr::pull(p_value_Fisher_E_J) %>%
        Rfast::bc()
    lambda_stat_all <- data %>% 
        dplyr::pull(stat_Fisher_all) %>%
        Rfast::bc()
    data %<>%
        dplyr::mutate(
            stat_Fisher_E_J_bc = ((stat_Fisher_E_J^lambda_stat_E_J)-1)/lambda_stat_E_J,
            stat_Fisher_all_bc = ((stat_Fisher_all^lambda_stat_all)-1)/lambda_stat_all,
            p_value_Fisher_bc = ((p_value_Fisher_E_J^lambda_p)-1)/lambda_p, 
            p_value_Fisher_lg = log(p_value_Fisher_E_J)) %>%
        dplyr::rename(p_value_Fisher = p_value_Fisher_E_J) %>% 
        dplyr::select(-c(p_value_Fisher_all))
    assign("lambda_p", lambda_p, envir = .GlobalEnv)
    assign(paste0("lambda_stat_E_J_", case), 
           lambda_stat_E_J, envir = .GlobalEnv)
    assign(paste0("lambda_stat_all_", case), 
           lambda_stat_all, envir = .GlobalEnv)
    return(data)
}

# Create tables for models
table_fun <- function(data_case, test.type) {
    call <- paste0("calls_", test.type)
    expand_grid(calls = get(call), expo) %>%
        # functional call, merge of power and call
        dplyr::mutate(dplyr::across(calls, str_replace_all, "power", 
                                    as.character(.$expo), .names = "formula")) %>%
        dplyr::mutate(purrr::map_dfr(formula, 
                                     function(x) RcppEigen::fastLm(formula(x), data = data_case_1) %>% lm_eval(data_case_1)))
}

Rcpp::cppFunction('
double RMSE_c (NumericVector pred, NumericVector dep) {
    double value = sqrt(sum(pow(pred - dep, 2))/pred.size());
    return value;
                  }')

# function automated lm    
lm_eval <- function(mod, data) {
    dep_var <- mod$formula[[2]]
    fitted_values <- mod$fitted.values
    mod$residuals <- c()
    mod$fitted.values <- c()
    
    values <- tibble(
        PRED = if(stringr::str_detect(dep_var, '_bc')){ 
            invBoxCox(fitted_values)
        } else if(stringr::str_detect(dep_var, '_lg')){
            exp(fitted_values)
        } else {fitted_values}, # y_hat 
        dependent = data %>% dplyr::pull(dep_var)
    )
    
    # computing corrected predictions
    values %<>% dplyr::mutate(PRED_cor = case_when(
        PRED <= 0 ~ 1e-12,
        PRED >= 1 ~ 1 - 1e-12,
        TRUE ~ PRED))
    
    # RMSE and corrected RMSE on full dataset
    RMSE <- RMSE_c(values$PRED, values$dependent)
    RMSE_cor <- RMSE_c(values$PRED_cor, values$dependent)
    
    # subset dataset: only the interesting part
    values_0.2 <- values %>%
        dplyr::filter(dependent >= 0.2)
    RMSE_0.2 <- RMSE_c(values_0.2$PRED, values_0.2$dependent)
    RMSE_cor_0.2 <- RMSE_c(values_0.2$PRED_cor, values_0.2$dependent)
    
    tibble(model  = list(mod),
           RMSE = RMSE,
           RMSE_cor = RMSE_cor,
           RMSE_0.2 = RMSE_0.2, 
           RMSE_cor_0.2 = RMSE_cor_0.2
    )
}

