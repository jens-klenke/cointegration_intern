#---- Preliminary ---- 
# packages 
source(here::here('01_code/packages/packages.R'))

# functions 
source(here::here('09_simulation_and_approximation-cdf/func_p_approximation.R'))

# load models 
load(here::here('09_simulation_and_approximation-cdf/models_package.RData'))

# lambda for box-cox
load(here::here('09_simulation_and_approximation-cdf/lambda_package.RData'))

# vom test selbst geliefert 
#k = 2 
#bh.test <- 21.1
#trendtype  <- 1
#test.type <- 'all'

get_p_value <- function(bh.test, k, trendtype, test.type){
    
    new_data <- tibble(stat_Fisher_all = bh.test, 
                   stat_Fisher_E_J = bh.test, 
                   k = k) %>%
        dplyr::mutate(stat_Fisher_E_J_bc = ((stat_Fisher_E_J^get_lambda(trendtype, 'stat', 'e_j'))-1)/get_lambda(trendtype, 'stat', 'e_j'),
                  stat_Fisher_all_bc = ((stat_Fisher_all^get_lambda(trendtype, 'stat', 'all'))-1)/get_lambda(trendtype, 'stat', 'all'))
    p.value_raw <- suppressWarnings(predict(get_model(trendtype, test.type), new_data))

    p.value_trans <- if(get_p_trans(trendtype, test.type) == 'log') {
        exp(p.value_raw)
        } else {
            if(get_p_trans(trendtype, test.type) == 'bc'){
                invBoxCox(p.value_raw)
            } else {p.value_raw}
        }

    p.value <- ifelse(p.value_trans >= 1, 9.9999e-1, ifelse(p.value_trans <= 0, 1e-12, p.value_trans))
    return(p.value)
}














#---- Test Plots ----
p_value_test_fun <- function(test.type){
    p_value_fun <- function(test_stat = 1:1000, trendtype, test.type, k){
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

