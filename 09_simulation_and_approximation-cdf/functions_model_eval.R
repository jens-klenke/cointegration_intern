### functions
best_5_table <- function(data, col){
    formattable(data%>%
                    dplyr::slice_min(RMSE_cor_0.2, n = 5)%>%
                    dplyr::select(-c(all_of(col), expo, model)),
                list(
                    RMSE = color_tile("green", "red"),
                    RMSE_cor = color_tile("green", "red"),
                    RMSE_0.2 = color_tile("green", "red"),
                    RMSE_cor_0.2 = color_tile("green", "red")
                ))
}

best_model <- function(data){
    data %<>% 
        dplyr::filter(RMSE_cor_0.2 == min(RMSE_cor_0.2)) %>%
        dplyr::select(formula)
}

invBoxCox <- function(x){
    x <- if (lambda_p == 0) exp(as.complex(x)) else (lambda_p*as.complex(x) + 1)^(1/lambda_p)
    return(Re(x))
}

# getting lambda
get_lambda <- function(data, case_w, art){
    data %>%
        dplyr::filter(case == case_w, 
                      side == art) %>%
        dplyr::select(value) %>%
        dplyr::pull()
}

# getting model 
get_model <- function(data, case_w, art){
    data %>%
        dplyr::filter(case == case_w,
                      test == art) %>%
        dplyr::select(models) %>%
        dplyr::pull() %>%
        purrr::pluck(1)
}

# getting data ready
# pred function must be changed  
plot_data <- function(data, case_w, art){
    
    data %>%
        dplyr::filter(case == case_w) %>%
        dplyr::mutate(stat_Fisher_E_J_bc = ((stat_Fisher_E_J^get_lambda(lambda_bc_EJ, case_w, 'stat'))-1)/get_lambda(lambda_bc_EJ, case_w, 'stat'),
                      stat_Fisher_all_bc = ((stat_Fisher_all^get_lambda(lambda_bc_ALL, case_w, 'stat'))-1)/get_lambda(lambda_bc_ALL, case_w, 'stat')) %>%
        modelr::add_predictions(get_model(final_models, case_w = case_w, art = art)) %>%
        dplyr::mutate(PRED = if(get_p_trans(final_models, case_w, art) == 'log'){exp(pred)} else{ invBoxCox(pred)}) %>%
        dplyr::mutate(PRED_cor = case_when(
            PRED <= 0 ~ 1e-12,
            PRED >= 1 ~ 1 - 1e-12,
            TRUE ~ PRED))
}

# getting transformation of the response
get_p_trans <- function(data, case_w, art){
    data %>%
        dplyr::filter(test == art,
                      case == case_w) %>%
        dplyr::pull(response)
}


# plot with facet k 
own_plot <- function(data, max_graph = 1){
    data %>%
        ggplot(aes(x = p_value_Fisher, y = PRED_cor)) +
        geom_line(color = '#004c93') +
        geom_segment(aes(x = 0, xend = max_graph, y = 0, yend = max_graph), 
                     linetype = 'dashed', size = 1, color = 'grey') +
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

