# packages 
source(here::here('01_code/packages/packages.R'))

# loading data 
load(file = here::here('09_simulation_and_approximation-cdf/plot_data_correctur.RData'))

# functions 
# plot with facet k 
own_plot <- function(data, max_graph = 1){
    data %>%
        dplyr::mutate(k = factor(k, labels = c('K = 1', 'K = 2', 'K = 3', 'K = 4', 'K = 5', 'K = 6',
                                               'K = 7', 'K = 8', 'K = 9', 'K = 10', 'K = 11')),
                      case = factor(trendtype, labels = c('Case = 1', 'Case = 2', 'Case = 3'))) %>%
        ggplot(aes(x = p_value_Fisher, y = PRED)) +
        geom_line(color = '#004c93') +
        xlim(c(0, max_graph))+
        ylim(c(0, max_graph))+
        labs(x = '\n Simulated p-values', y = 'Approximated p-values \n')+
        theme_bw()+
        facet_grid(k ~ case)+
        scale_x_continuous(sec.axis = sec_axis(~.*10)) +
        scale_y_continuous(sec.axis = sec_axis(~.*10), 
                           breaks = seq(0, 1, 0.5)) +
        theme(panel.spacing = unit(1, "lines"),
              strip.background = element_rect(colour = 'black',
                                              fill = '#004c93'),
              strip.text.x = element_text(color = 'white'), 
              strip.text.y = element_text(color = 'white'), 
              axis.ticks.y.right = element_line(colour = 'white'), 
              axis.text.y.right =  element_text(colour = "white"),
              axis.ticks.x.top = element_line(colour = 'white'), 
              axis.text.x.top =  element_text(colour = "white")
        )
}

own_plot_0.2 <- function(data, max_graph = 0.2){
    data %>%
        dplyr::mutate(k = factor(k, labels = c('K = 1', 'K = 2', 'K = 3', 'K = 4', 'K = 5', 'K = 6',
                                               'K = 7', 'K = 8', 'K = 9', 'K = 10', 'K = 11')),
                      case = factor(trendtype, labels = c('Case = 1', 'Case = 2', 'Case = 3'))) %>%
        dplyr::filter(p_value_Fisher <= 0.2) %>%
        ggplot(aes(x = p_value_Fisher, y = PRED)) +
        geom_line(color = '#004c93') +
        xlim(c(0, max_graph))+
        ylim(c(0, max_graph))+
        labs(x = '\n Simulated p-values', y = 'Approximated p-values \n')+
        theme_bw()+
        facet_grid(k ~ case) +
        scale_x_continuous(sec.axis = sec_axis(~.*10),
                           breaks = seq(0, 1, 0.1)) +
        scale_y_continuous(sec.axis = sec_axis(~.*10), 
                           breaks = seq(0, 1, 0.1)) +
        theme(panel.spacing = unit(1, "lines"),
              strip.background = element_rect(colour = 'black',
                                              fill = '#004c93'),
              strip.text.x = element_text(color = 'white'), 
              strip.text.y = element_text(color = 'white'), 
              axis.ticks.y.right = element_line(colour = 'white'), 
              axis.text.y.right =  element_text(colour = "white"),
              axis.ticks.x.top = element_line(colour = 'white'), 
              axis.text.x.top =  element_text(colour = "white")
        )
}




#-- plot sim-p-values vs predicted p-values ----

# diagnose plots
p.sim_p.aprox_all <- own_plot(data_all)
p.sim_p.aprox_all_0.2 <- own_plot_0.2(data_all)

# diagnose plots
p.sim_p.aprox_e_j <- own_plot(data_e_j)
p.sim_p.aprox_e_j_0.2 <- own_plot_0.2(data_e_j)

# saving plots 
save(p.sim_p.aprox_all, p.sim_p.aprox_all_0.2, p.sim_p.aprox_e_j, p.sim_p.aprox_e_j_0.2,
     plot_p_stat_e_j, plot_p_stat_all,
     file = here::here('09_simulation_and_approximation-cdf/01_paper/paper_plots.RData'))