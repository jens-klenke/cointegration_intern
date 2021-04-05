# functions

# function automated lm    
own_lm <- function(call_mod, data){
    lm(paste(call_mod, collapse = ''),  data = data)
}

# inverse BoxCox function
invBoxCox <- function(x){
    x <- if (lambda_p == 0) exp(as.complex(x)) else (lambda_p*as.complex(x) + 1)^(1/lambda_p)
    return(Re(x))
}

# BoxCox + log Transformation
bc_log_E_J_fun <- function(data, lambda_stat) {
    data %>%
        mutate(
            stat_Fisher_E_J_bc = ((stat_Fisher_E_J^lambda_stat)-1)/lambda_stat,
            p_value_Fisher_E_J_bc = ((p_value_Fisher_E_J^lambda_p)-1)/lambda_p, 
            p_value_Fisher_E_J_lg = log(p_value_Fisher_E_J))
}

# metric function
metric_fun <- function(object, data){
    
    # save call and case
 #   mod_call <- paste(colnames(object$model)[1], '~', paste(colnames(object$model)[2:length(colnames(object$model))], collapse =' '))
    
#    case <- object$call$data
    
    # dependent variable 
    dep_var <- colnames(object$model)[1]
    
    # dataset for prediction
    values <- tibble(
        PRED = if(str_detect(dep_var, '_bc')){ 
            invBoxCox(object$fitted.values)
        } else if(str_detect(dep_var, '_lg')){
            exp(object$fitted.values)
        } else {object$fitted.values}, # y_hat 
        dependent = if(str_detect(dep_var, 'E_J')){ 
            data$p_value_Fisher_E_J
        } else {data$p_value_Fisher_all},
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
        dplyr::filter(dependent >= 0.8)
    RMSE_0.2 <- sqrt((sum((values_0.2$PRED - values_0.2$dependent)^2)/nrow(values_0.2)))
    RMSE_cor_0.2 <- sqrt((sum((values_0.2$PRED_cor - values_0.2$dependent)^2)/nrow(values_0.2)))
    
    # values for the plots 
    values <- values%>%
        dplyr::filter(dependent %in% seq(0.001, 1, 0.001))
    
    mod_sum <- tibble(#model = as.character(mod_call),
                      RMSE = as.numeric(RMSE),
                      RMSE_cor = as.numeric(RMSE_cor),
                      RMSE_0.2 = as.numeric(RMSE_0.2), 
                      RMSE_cor_0.2 = as.numeric(RMSE_cor_0.2), 
#                      case = as.character(case),
                      pred = list(values)
    )
    
    return(mod_sum)
}

# Create tables for models
table_E_J_fun <- function(data_case) {
    expand_grid(calls_E_J, expo) %>%
        # functional call, merge of power and call
        dplyr::mutate(dplyr::across("calls_E_J", str_replace_all, "power", 
                                    as.character(.$expo), .names = "formula")) %>%
        # adding data
        dplyr::mutate(data = list(data_case)) %>%
        # fitting the model 
        dplyr::mutate(models = map2(formula, data, own_lm)) %>%
        # calculating the metrics for model evaluation
        dplyr::mutate(map2_df(models, data, metric_fun)) %>%
        # deleting data and other unimportant variables
        dplyr::select(-data)
}

table_all_fun <- function(data_case) {
    expand_grid(calls_all, expo) %>%
        # functional call, merge of power and call
        dplyr::mutate(dplyr::across("calls_all", str_replace_all, "power", 
                                    as.character(.$expo), .names = "formula")) %>%
        # adding data
        dplyr::mutate(data = list(data_case)) %>%
        # fitting the model 
        dplyr::mutate(models = map2(formula, data, own_lm)) %>%
        # calculating the metrics for model evaluation
        dplyr::mutate(map2_df(models, data, metric_fun)) %>%
        # deleting data and other unimportant variables
        dplyr::select(-data)
}
