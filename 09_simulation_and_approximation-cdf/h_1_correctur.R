# packages 
source(here::here('01_code/packages/packages.R'))

# functions 
source(here::here('09_simulation_and_approximation-cdf/func_p_approximation.R'))
source(here::here('09_simulation_and_approximation-cdf/functions_package.R'))

# load model matrix
load(here::here('09_simulation_and_approximation-cdf/models.RData'))

# H_1 limit 
load(here::here('09_simulation_and_approximation-cdf/H_1-values.RData'))

# parallel 
plan(multisession, workers = 3)

# functions 
critical_fun <- function(case_s, test_s){
    crit_val %>%
        dplyr::filter(case == case_s, 
                      test == test_s) %>%
        dplyr::ungroup() %>%
        dplyr::select(-c(case, test))
}

critical_fun_sim <- function(case_s, test_s){
    crit_val %>%
        dplyr::filter(case == case_s, 
                      test == test_s) %>%
        dplyr::ungroup() %>%
        dplyr::select(-c(case, test))
}

# storing all critical_val in one tibble 
crit_val <- dplyr::bind_rows(
    # for the all test
    crit_val_all_100 %>%
        dplyr::rename(crit_val = 'stat_Fisher_all') %>%
        dplyr::mutate(test = 'all'),
    
    # for the e_j test
    crit_val_e_j_100 %>%
        dplyr::mutate(test = 'E_J') %>%
        dplyr::rename(crit_val = 'stat_Fisher_E_J'))

#-- graphical inspection ----

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

# predicting the p-values
p_value_plots_e_j <- p_value_test_fun('E_J')
p_value_plots_all <- p_value_test_fun('all')


# getting minimum bh.test value bei dem p.vlaue wieder steigt

h_1_cor_min <- bind_rows( p_value_plots_all %>%
    dplyr::group_by(k, trendtype) %>%
    dplyr::filter(p_value == min(p_value)) %>%
    dplyr::filter(test_stat == min(test_stat)) ,

p_value_plots_e_j %>% 
    dplyr::group_by(k, trendtype) %>%
    dplyr::filter(p_value == min(p_value)) %>%
    dplyr::filter(test_stat == min(test_stat)))


critical_fun_sim <- function(case_s, test_s){
    h_1_cor_min %>%
        dplyr::rename(crit_val = test_stat) %>%
        dplyr::filter(trendtype == case_s, 
                      test.type == test_s) %>%
        dplyr::ungroup() %>%
        dplyr::select(c(k, crit_val))
}


# adding critical val to tibble for last 100 value 
#models %<>%
#    dplyr::mutate(critical = purrr::map2(case, test.type, critical_fun))

# adding critical val to tibble
models %<>%
    dplyr::mutate(critical = purrr::map2(case, test.type, critical_fun_sim))

# save model matrix with critcal values
# save(models, file = here::here('09_simulation_and_approximation-cdf/models.RData'))

test_plot_fun(p_value_plots_e_j)
test_plot_fun(p_value_plots_all)

#-- plot sim-p-values vs predicted p-values ----

# Load Data
#if(Sys.info()['nodename'] != "DELL-ARBEIT") {
#    load(here::here('09_simulation_and_approximation-cdf/data_cases.RData'))
#} else if(Sys.info()['nodename'] == "DELL-ARBEIT") { # Jens 
#    load('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\data_cases.RData')
#}

load(file = here::here('09_simulation_and_approximation-cdf/plot_data_correctur.RData'))

#data <- bind_rows(data_case_1, data_case_2, data_case_3)

#rm(data_case_1, data_case_2, data_case_3)

#data %<>%
#    dplyr::filter(p_value_Fisher %in% c(min(p_value_Fisher), seq(0, 1, 0.0001)))

#data_all <- data %>%
#    dplyr::rename(bh.test = stat_Fisher_all,
#                  trendtype = case) %>%
#    dplyr::mutate(test.type = 'all') %>%
#    dplyr::select(bh.test, trendtype, test.type, k, p_value_Fisher) %>%
#    dplyr::mutate(PRED = furrr::future_pmap_dbl(., get_p_value))

#data_e_j <- data %>%
#    dplyr::rename(bh.test = stat_Fisher_E_J,
#                  trendtype = case) %>%
#    dplyr::mutate(test.type = 'E_J') %>%
#    dplyr::select(bh.test, trendtype, test.type, k, p_value_Fisher) %>%
#    dplyr::mutate(PRED = furrr::future_pmap_dbl(., get_p_value))


# plot with facet k 
own_plot <- function(data, max_graph = 1){
    data %>%
        ggplot(aes(x = p_value_Fisher, y = PRED)) +
        geom_line(color = '#004c93') +
        xlim(c(0, max_graph))+
        ylim(c(0, max_graph))+
        labs(x = '\n Simulated p-values', y = 'Approximated p-values \n')+
        theme_bw()+
        facet_grid(k ~ trendtype)+
        scale_y_continuous(sec.axis = sec_axis(~.*10, name = "k"), 
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
        dplyr::filter(p_value_Fisher <= 0.2) %>%
        ggplot(aes(x = p_value_Fisher, y = PRED)) +
        geom_line(color = '#004c93') +
        xlim(c(0, max_graph))+
        ylim(c(0, max_graph))+
        labs(x = '\n Simulated p-values', y = 'Approximated p-values \n')+
        theme_bw()+
        facet_grid(k ~ trendtype)+
        scale_y_continuous(sec.axis = sec_axis(~.*10, name = "k"), 
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



# diagnose plots
own_plot(data_all)
own_plot_0.2(data_all)

# diagnose plots
own_plot(data_e_j)
own_plot_0.2(data_e_j)

# save(data_all, data_e_j, file = here::here('09_simulation_and_approximation-cdf/plot_data_correctur.RData'))




# plot
p_value_plots_e_j %>%
    tidyr::pivot_longer(cols = c(p_value, p_value_2), names_to =  c('variable')) %>%
    ggplot2::ggplot(aes(x = test_stat, y = value, colour = variable)) +
    geom_line(size = 1) +
    theme_bw() +
    facet_grid(k ~ trendtype) +
    scale_fill_manual(values = c('#004c93', '#f51137')) +
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










