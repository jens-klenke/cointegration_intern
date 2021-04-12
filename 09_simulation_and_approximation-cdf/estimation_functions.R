# functions

# function automated lm    
own_lm <- function(call_mod, data){
    # lm fit
    mod <- lm(call_mod,  data = data)
    # dependent variable 
    dep_var <- colnames(mod$model)[1]
    # fitted values
    fitted_values <- mod$fitted.values
    # cleand model
    mod <- clean_lm(mod)
    
    tib <-tibble(model  = list(mod), 
                 dep_var = dep_var,
                 fitted_values = list(fitted_values))
    
    return(tib)
}

# clean model 
clean_lm <- function(object) {
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
    
    object$family$variance = c()
    object$family$dev.resids = c()
    object$family$aic = c()
    object$family$validmu = c()
    object$family$simulate = c()
    attr(object$terms,".Environment") = c()
    attr(object$formula,".Environment") = c()
    
    object
}

# inverse BoxCox function
invBoxCox <- function(x){
    x <- if (lambda_p == 0) exp(as.complex(x)) else (lambda_p*as.complex(x) + 1)^(1/lambda_p)
    return(Re(x))
}

# BoxCox + log Transformation
bc_log_fun <- function(data) {
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
    return(data)
}



# Create tables for models
table_E_J_fun <- function(data_case) {
    data <- data_case
    expand_grid(calls_E_J, expo) %>%
        # functional call, merge of power and call
        dplyr::mutate(dplyr::across("calls_E_J", str_replace_all, "power", 
                                    as.character(.$expo), .names = "formula")) %>%
        # fitting the model 
        dplyr::mutate(map_df(formula, ~own_lm(call_mod = ., data = data))) %>%
        # calculating the metrics for model evaluation
        dplyr::mutate(furrr::future_pmap_dfr(list(fitted_values, dep_var), 
                              ~new_metric_fun(.x, .y, data))) %>%
        # deleting data and other unimportant variables
        dplyr::select(-fitted_values)
}

table_all_fun <- function(data_case) {
    data <- data_case
    expand_grid(calls_all, expo) %>%
        # functional call, merge of power and call
        dplyr::mutate(dplyr::across("calls_all", str_replace_all, "power", 
                                    as.character(.$expo), .names = "formula")) %>%
        # fitting the model 
        dplyr::mutate(map_df(formula, ~own_lm(call_mod = ., data = data))) %>%
        # calculating the metrics for model evaluation
        dplyr::mutate(furrr::future_pmap_dfr(list(fitted_values, dep_var), 
                              ~new_metric_fun(.x, .y, data))) %>%
        # deleting data and other unimportant variables
        dplyr::select(-fitted_values)
}



# metric function
new_metric_fun <- function(fitted_values, dep_var, data, ...){
    # dataset for prediction
    values <- tibble(
        PRED = if(str_detect(dep_var, '_bc')){ 
            invBoxCox(fitted_values)
        } else if(str_detect(dep_var, '_lg')){
            exp(fitted_values)
        } else {fitted_values}, # y_hat 
        dependent = data$p_value_Fisher,
        k = data$k
    )
    
    
    # computing corrected predictions
    values <- values%>%
        dplyr::mutate(PRED_cor = case_when(
            PRED <= 0 ~ 1e-12,
            PRED >= 1 ~ 1 - 1e-12,
            TRUE ~ PRED))
    
    # RMSE and corrected RMSE on the full dataset
    RMSE <- sqrt((sum((values$PRED - values$dependent)^2)/nrow(values)))
    RMSE_cor <- sqrt((sum((values$PRED_cor - values$dependent)^2)/nrow(values)))
    
    # smaller dataset only the intersting part, correction is at 
    values_0.2 <- values %>%
        dplyr::filter(dependent >= 0.2)
    RMSE_0.2 <- sqrt((sum((values_0.2$PRED - values_0.2$dependent)^2)/nrow(values_0.2)))
    RMSE_cor_0.2 <- sqrt((sum((values_0.2$PRED_cor - values_0.2$dependent)^2)/nrow(values_0.2)))
    
    mod_sum <- tibble(#model = as.character(mod_call),
        RMSE = as.numeric(RMSE),
        RMSE_cor = as.numeric(RMSE_cor),
        RMSE_0.2 = as.numeric(RMSE_0.2), 
        RMSE_cor_0.2 = as.numeric(RMSE_cor_0.2) 
    )
    
    return(mod_sum)
}

