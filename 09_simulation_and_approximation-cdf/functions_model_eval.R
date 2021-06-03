### functions
best_5_table <- function(data){
    formattable(data%>%
                    dplyr::slice_min(RMSE_cor_0.2, n = 5)%>%
                    dplyr::select(-c(calls, expo, model)),
                list(
                    RMSE = color_tile("green", "red"),
                    RMSE_cor = color_tile("green", "red"),
                    RMSE_0.2 = color_tile("green", "red"),
                    RMSE_cor_0.2 = color_tile("green", "red")
                ))
}

invBoxCox <- function(x){
    x <- if (lambda_p == 0) exp(as.complex(x)) else (lambda_p*as.complex(x) + 1)^(1/lambda_p)
    return(Re(x))
}

# getting lambda
get_lambda <- function(data, case_w, art, test_w){
    data %>%
        dplyr::filter(case == case_w, 
                      side == art, 
                      test == test_w) %>%
        dplyr::select(value) %>%
        dplyr::pull()
}

# getting data ready
# pred function must be changed  
add_pred <- function(data, art){
    
    # extract case
    case <- deparse(substitute(data)) %>% stringr::str_extract("[0-9]+")
    
    # extract model 
    model <- get_model_eval(get(paste0('table_', art, '_case_', case))) %>%
        add_predvars(get(paste0("data_case_", case)))
    
    data %>%
        dplyr::mutate(PRED =  as.vector(model.matrix(model$terms, data) %*% coef(model))) %>%
        dplyr::mutate(PRED = if(stringr::str_detect(get_p_trans_eval(model), '_lg')){
            exp(PRED)
        } else if(stringr::str_detect(get_p_trans_eval(model), '_bc')){
            invBoxCox(PRED)
        } else {
            PRED
        }) %>%
        dplyr::mutate(PRED_cor = case_when(
            PRED <= 0 ~ 1e-12,
            PRED >= 1 ~ 1 - 1e-12,
            TRUE ~ PRED))
}


get_p_trans_eval <- function(model){
    model %>%
        purrr::pluck('formula') %>%
        purrr::pluck(2)
}

get_model_eval <- function(data){
    data %>% 
        dplyr::slice_min(RMSE_0.2, n = 1) %>%
        dplyr::select(model) %>%
        dplyr::pull() %>%
        purrr::pluck(1)
}

# plot with facet k 
own_plot <- function(data, max_graph = 1){
    data %>%
        ggplot(aes(x = p_value_Fisher, y = PRED_cor)) +
        geom_line(color = '#004c93') +
        xlim(c(0, max_graph))+
        ylim(c(0, max_graph))+
        labs(x = '\n Simulated p-values', y = 'Approximated p-values \n')+
        theme_bw()+
        facet_wrap(~k)+
        theme(panel.spacing = unit(1, "lines"),
              strip.background = element_rect(colour = 'black',
                                              fill = '#004c93'),
              strip.text.x = element_text(size = 12, color = 'white' # , face = "bold.italic"
              ),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20))
}

own_plot_0.2 <- function(data, max_graph = 0.2){
    data %>%
        dplyr::filter(p_value_Fisher <= 0.2) %>%
        ggplot(aes(x = p_value_Fisher, y = PRED_cor)) +
        geom_line(color = '#004c93') +
        xlim(c(0, max_graph))+
        ylim(c(0, max_graph))+
        labs(x = '\n Simulated p-values', y = 'Approximated p-values \n')+
        theme_bw()+
        facet_wrap(~k)+
        theme(panel.spacing = unit(1, "lines"),
              strip.background = element_rect(colour = 'black',
                                              fill = '#004c93'),
              strip.text.x = element_text(size = 12, color = 'white' # , face = "bold.italic"
              ),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20))
}

# add terms to models
add_predvars <- function(object, data) {
    object$terms <- terms(object$formula)
    vars <- attr(object$terms, "variables")
    variables <- eval(vars, data)
    varnames <- vars %>% length() - 1
    predvars <- vars
    
    for (i in 1:varnames) predvars[[i + 1L]]  <- makepredictcall(variables[[i]], vars[[i + 1L]])
    attr(object$terms, "predvars") <- predvars
    object
}

best_model_fun <- function(art, case){
    # extract best model and add predvars
    get_model_eval(get(paste0('table_', art, '_case_', case))) %>%
        add_predvars(get(paste0("data_case_", case)))
}