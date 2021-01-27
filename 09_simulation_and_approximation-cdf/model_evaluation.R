#-- functions ----
pred_1000 <- function(data, object, dep_var){
    model <- deparse(substitute(object))
    
    data <- data%>%
        dplyr::filter(!!rlang::sym(dep_var) %in% seq(0.000, 1, 0.001))%>%
        do(modelr::add_predictions(., object))%>%
        dplyr::select(pred)%>%
        dplyr::rename(!!paste0('PRED_', model) := pred)
    
    return(data)
}

k_com <- function(data){
    
    dep <- data [,'dependent']%>%
        pull()
    K = 1
    
    k <- rep(NA, length(dep)) 
    k[1] <- 1
    
    for (i in 2:length(dep)){
        
        if(((dep[i-1] - dep[i]) < 0) == FALSE){
            K = K + 1} 
        
        k[i] <- K
    } 
    return(k)
} 

#---- Preliminary ---- 
source(here::here('01_code/packages/packages.R'))

# load metrics
load(here::here("09_simulation_and_approximation-cdf/model_metrics_ALL_server.Rdata"))

load(here::here("09_simulation_and_approximation-cdf/model_metrics_E_J_server.Rdata"))

#-- Analysis ALL ----

# tables 

table_all_case_1 <- formattable(model_metrics_ALL%>%
                                    dplyr::filter(case == 'data_case_1')%>%
                                    dplyr::slice_min(RMSE_0.2, n = 5)%>%
                                    dplyr::select(-c(pred, case)),
                                list(
                                    RMSE = color_tile("green", "red"),
                                    RMSE_cor = color_tile("green", "red"),
                                    RMSE_0.2 = color_tile("green", "red"),
                                    RMSE_cor_0.2 = color_tile("green", "red")
))

table_all_case_2 <- formattable(model_metrics_ALL%>%
                                    dplyr::filter(case == 'data_case_2')%>%
                                    dplyr::slice_min(RMSE_0.2, n = 5)%>%
                                    dplyr::select(-c(pred, case)),
                                list(
                                    RMSE = color_tile("green", "red"),
                                    RMSE_cor = color_tile("green", "red"),
                                    RMSE_0.2 = color_tile("green", "red"),
                                    RMSE_cor_0.2 = color_tile("green", "red")
                                ))

table_all_case_3 <- formattable(model_metrics_ALL%>%
                                    dplyr::filter(case == 'data_case_3')%>%
                                    dplyr::slice_min(RMSE_0.2, n = 5)%>%
                                    dplyr::select(-c(pred, case)),
                                list(
                                    RMSE = color_tile("green", "red"),
                                    RMSE_cor = color_tile("green", "red"),
                                    RMSE_0.2 = color_tile("green", "red"),
                                    RMSE_cor_0.2 = color_tile("green", "red")
                                ))


# plots 
plot_data_ALL_case_1 <- 
        model_metrics_ALL%>%
    dplyr::filter(case == 'data_case_1')%>%
    dplyr::filter(RMSE == min(RMSE))%>%
    unnest(pred)%>%
    dplyr::mutate(k = k_com(.))%>%
    dplyr::mutate(k = factor(.$k, levels = c(1:11), 
                             labels = c(paste('k =', 1:11))))%>%
    ggplot(aes(x = dependent, y = PRED_cor))+
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

plot_data_ALL_case_2 <- 
    model_metrics_ALL%>%
    dplyr::filter(case == 'data_case_2')%>%
    dplyr::filter(RMSE == min(RMSE))%>%
    unnest(pred)%>%
    dplyr::mutate(k = k_com(.))%>%
    dplyr::mutate(k = factor(.$k, levels = c(1:11), 
                             labels = c(paste('k =', 1:11))))%>%
    ggplot(aes(x = dependent, y = PRED_cor))+
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

plot_data_ALL_case_3 <- 
    model_metrics_ALL%>%
    dplyr::filter(case == 'data_case_3')%>%
    dplyr::filter(RMSE == min(RMSE))%>% # condition may need to be changes
    unnest(pred)%>%
    dplyr::mutate(k = k_com(.))%>%
    dplyr::mutate(k = factor(.$k, levels = c(1:11), 
                             labels = c(paste('k =', 1:11))))%>%
    ggplot(aes(x = dependent, y = PRED_cor))+
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

#-- Analysis E_J ----

# tables 

table_e_j_case_1 <- formattable(model_metrics_E_J%>%
                                    dplyr::filter(case == 'data_case_1')%>%
                                    dplyr::slice_min(RMSE_0.2, n = 5)%>%
                                    dplyr::select(-c(pred, case)),
                                list(
                                    RMSE = color_tile("green", "red"),
                                    RMSE_cor = color_tile("green", "red"),
                                    RMSE_0.2 = color_tile("green", "red"),
                                    RMSE_cor_0.2 = color_tile("green", "red")
                                ))

table_e_j_case_2 <- formattable(model_metrics_E_J%>%
                                    dplyr::filter(case == 'data_case_2')%>%
                                    dplyr::slice_min(RMSE_0.2, n = 5)%>%
                                    dplyr::select(-c(pred, case)),
                                list(
                                    RMSE = color_tile("green", "red"),
                                    RMSE_cor = color_tile("green", "red"),
                                    RMSE_0.2 = color_tile("green", "red"),
                                    RMSE_cor_0.2 = color_tile("green", "red")
                                ))

table_e_j_case_3 <- formattable(model_metrics_E_J%>%
                                    dplyr::filter(case == 'data_case_3')%>%
                                    dplyr::slice_min(RMSE_0.2, n = 5)%>%
                                    dplyr::select(-c(pred, case)),
                                list(
                                    RMSE = color_tile("green", "red"),
                                    RMSE_cor = color_tile("green", "red"),
                                    RMSE_0.2 = color_tile("green", "red"),
                                    RMSE_cor_0.2 = color_tile("green", "red")
                                ))


# plots 
plot_data_e_j_case_1 <- 
    model_metrics_E_J%>%
    dplyr::filter(case == 'data_case_1')%>%
    dplyr::filter(RMSE == min(RMSE))%>%
    unnest(pred)%>%
    dplyr::mutate(k = k_com(.))%>%
    dplyr::mutate(k = factor(.$k, levels = c(1:11), 
                             labels = c(paste('k =', 1:11))))%>%
    ggplot(aes(x = dependent, y = PRED_cor))+
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

plot_data_e_j_case_2 <- 
    model_metrics_E_J%>%
    dplyr::filter(case == 'data_case_2')%>%
    dplyr::filter(RMSE == min(RMSE))%>%
    unnest(pred)%>%
    dplyr::mutate(k = k_com(.))%>%
    dplyr::mutate(k = factor(.$k, levels = c(1:11), 
                             labels = c(paste('k =', 1:11))))%>%
    ggplot(aes(x = dependent, y = PRED_cor))+
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

plot_data_e_j_case_3 <- 
    model_metrics_E_J%>%
    dplyr::filter(case == 'data_case_3')%>%
    dplyr::filter(RMSE == min(RMSE))%>% # condition may need to be changes
    unnest(pred)%>%
    dplyr::mutate(k = k_com(.))%>%
    dplyr::mutate(k = factor(.$k, levels = c(1:11), 
                             labels = c(paste('k =', 1:11))))%>%
    ggplot(aes(x = dependent, y = PRED_cor))+
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
