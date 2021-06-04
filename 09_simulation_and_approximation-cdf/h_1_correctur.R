# packages 
source(here::here('01_code/packages/packages.R'))

# functions 
source(here::here('09_simulation_and_approximation-cdf/func_p_approximation.R'))

# load model matrix
load(here::here('09_simulation_and_approximation-cdf/models.RData'))

# H_1 limit 
load(here::here('09_simulation_and_approximation-cdf/H_1-values.RData'))

x# functions 
critical_fun <- function(case_s, test_s){
    crit_val %>%
        dplyr::filter(case == case_s, 
                      test == test_s) %>%
        dplyr::ungroup() %>%
        dplyr::select(-c(case, test))
}

# storing all critical_val in one tibble 
crit_val <- dplyr::bind_rows(
    # for the all test
    crit_val_all_10 %>%
        dplyr::rename(crit_val = 'stat_Fisher_all') %>%
        dplyr::mutate(test = 'all'),
    
    # for the e_j test
    crit_val_e_j_10 %>%
        dplyr::mutate(test = 'E_J') %>%
        dplyr::rename(crit_val = 'stat_Fisher_E_J'))


# adding critical val to tibble

models %<>%
    dplyr::mutate(critical = purrr::map2(case, test.type, critical_fun))

#-- graphical inspection ----




p_value_test_fun <- function(test.type){
    p_value_fun <- function(test_stat = 1:100, trendtype, test.type, k){
        tibble(
            test_stat = test_stat,
            p_value = purrr::map_dbl(test_stat, ~get_p_value(., trendtype, test.type, k)), 
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

test_plot_fun <- function(p_values){
    p_values %>%
        ggplot(aes(x = test_stat, y = p_value)) +
        geom_line(color = '#004c93') +
        labs(x = 'Test Statistic', y = 'Approximated p-values \n') +
        theme_bw() +
        facet_grid(k ~ trendtype) +
        scale_y_continuous(sec.axis = sec_axis(~.*10, name = "k"), 
                           breaks = seq(0, 1, 0.5)) +
        scale_x_continuous(sec.axis = sec_axis(~.*10, name = "Deterministic component")) +
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

p_value_plots_e_j <- p_value_test_fun('e_j')
p_value_plots_all <- p_value_test_fun('all')

test_plot_fun(p_value_plots_e_j)
test_plot_fun(p_value_plots_all)










p_values %>%
    ggplot(aes(x = test_stat, y = p_value)) +
    geom_line(color = '#004c93') +
    labs(x = 'Test Statistic', y = 'Approximated p-values \n') +
    