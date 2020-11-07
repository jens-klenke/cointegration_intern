##### function
pred_1000 <- function(data, object, dep_var){
    model <- deparse(substitute(object))
    
    data <- data%>%
        dplyr::filter(!!rlang::sym(dep_var) %in% seq(0.000, 1, 0.001))%>%
        do(modelr::add_predictions(., object))%>%
        dplyr::select(pred)%>%
        dplyr::rename(!!paste0('PRED_', model) := pred)
    
    return(data)
}


# data for the prediction
plot_data <- data_case_1%>%
    dplyr::filter(p_value_E_G %in% seq(0.000, 1, 0.001))%>%
    dplyr::mutate(k = factor(.$k, levels = c(1:11), 
                             labels = c(paste('k =', 1:11))))%>%
    # predictions
    dplyr::mutate(pred_1000(data_case_1, mod_E_G_case.1_p_4, 'p_value_E_G'),
                  pred_1000(data_case_1, mod_E_G_case.1_p_3, 'p_value_E_G'))

### ggplot 
plot_data%>%
#    dplyr::filter(k %in% 1:2)%>%
ggplot(aes(x = p_value_E_G, y = PRED_mod_E_G_case.1_p_4))+
    geom_line()+
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype = 'dashed', size = 1, color = 'grey')+
    xlim(c(0, 1))+
    ylim(c(0, 1))+
    labs(x = '\n Simulated p-values', y = 'Approximated p-values \n')+
    theme_bw()+
    facet_wrap(~k)+
    theme(panel.spacing = unit(1, "lines"), 
          strip.background = element_rect(colour="black",
                                          fill = "#004c93"), 
          strip.text.x = element_text(
              size = 12, color = "white" # , face = "bold.italic"
          ), 
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20))

     



 