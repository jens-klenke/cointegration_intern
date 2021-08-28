#-- plot functions ----
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

p_value_test_fun <- function(test.type){
    p_value_fun <- function(test_stat = seq(1, 100, 1), trendtype, test.type, k){
        tibble(
            test_stat = test_stat,
            p_value = purrr::map_dbl(test_stat, ~get_p_value(., trendtype, test.type, k)),
            p_value_2 = purrr::map_dbl(test_stat, ~get_p_value_2(., trendtype, test.type, k)),
            trendtype = trendtype, 
            test.type = test.type,
            k = k)
    }
    expand_grid(
        trendtype = 1:3, 
        test.type = test.type, 
        k = 1:11
    ) %>% 
        pmap_df(p_value_fun)
}

plot_p_stat <-function(data){
    data %>%
        tidyr::pivot_longer(cols = c(p_value, p_value_2), names_to =  c('variable')) %>%
        dplyr::mutate(k = factor(k, labels = c('K = 1', 'K = 2', 'K = 3', 'K = 4', 'K = 5', 'K = 6',
                                               'K = 7', 'K = 8', 'K = 9', 'K = 10', 'K = 11')),
                      case = factor(trendtype, labels = c('Case = 1', 'Case = 2', 'Case = 3'))) %>%
        ggplot2::ggplot(aes(x = test_stat, y = value, colour = variable)) +
        geom_line(size = 0.5) +
        theme_bw() + 
        labs(x = '\n Test Statistic', y = 'p-value \n') +
        scale_color_manual(values = c('#004c93', '#f51137'),
                           name = 'Approximations',
                           labels = c("Corrected", "Uncorrected")) +
        facet_grid(k ~ case) +
        scale_y_continuous(sec.axis = sec_axis(~.*10), 
                           breaks = seq(0, 1, 0.5),
                           limits = c(0, 1)) +
        scale_x_continuous(sec.axis = sec_axis(~.*10),
                           limits = c(0, 100)) +
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


plot_p_stat_k.1 <-function(data){
    data %>%
        tidyr::pivot_longer(cols = c(p_value, p_value_2), names_to =  c('variable')) %>%
        dplyr::mutate(k = factor(k, labels = c('K = 1', 'K = 2', 'K = 3', 'K = 4', 'K = 5', 'K = 6',
                                               'K = 7', 'K = 8', 'K = 9', 'K = 10', 'K = 11')),
                      case = factor(trendtype, labels = c('Case = 1', 'Case = 2', 'Case = 3'))) %>%
        dplyr::filter(k == 'K = 1') %>%
        ggplot2::ggplot(aes(x = test_stat, y = value, colour = variable)) +
        geom_line(size = 0.5) +
        theme_bw() + 
        labs(x = '\n Test Statistic', y = 'p-value \n') +
        scale_color_manual(values = c('#004c93', '#f51137'),
                           name = 'Approximations',
                           labels = c("Corrected", "Uncorrected")) +
        facet_grid(~ case) +
        scale_y_continuous(sec.axis = sec_axis(~.*10), 
                           breaks = seq(0, 1, 0.5),
                           limits = c(0, 1)) +
        scale_x_continuous(sec.axis = sec_axis(~.*10),
                           limits = c(0, 100)) +
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
