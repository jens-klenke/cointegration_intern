## Load Simulation Data 
load("/Users/Janine/Desktop/BayerHanck/data_cases.RData")

library(ggplot2)
library(dplyr)

dens_plot <- function(data, k, var){
    data %>% 
        dplyr::filter(k == k) %>%
        ggplot() + 
        geom_density(aes(x = get(var)), colour = '#004c93') +
        xlab("Test Statistic") +
        theme_bw() +
        theme(panel.spacing = unit(1, "lines"),
              strip.background = element_rect(colour = 'black',
                                              fill = '#004c93'),
              strip.text.x = element_text(color = 'white'), 
              strip.text.y = element_text(color = 'white'), 
              axis.ticks.y.right = element_line(colour = 'white'), 
              axis.text.y.right =  element_text(colour = "white"),
              axis.ticks.x.top = element_line(colour = 'white'), 
              axis.text.x.top =  element_text(colour = "white"),
              plot.title = element_text(hjust = 0.5)
              )
}

# all tests
ts_dens <- dens_plot(data_case_1, 1, "stat_Fisher_all")
ts_dens_bc <- dens_plot(data_case_1, 1, "stat_Fisher_all_bc")
dense_all <- gridExtra::grid.arrange(ts_dens + 
                                         ggtitle("No transformation"), 
                                     ts_dens_bc + ggtitle("Box-Cox transformation"), 
                                     ncol = 2)

# EJ
ts_dens_EJ <- dens_plot(data_case_1, 1, "stat_Fisher_E_J")
ts_dens_bc_EJ <- dens_plot(data_case_1, 1, "stat_Fisher_E_J_bc")
dense_EJ <- gridExtra::grid.arrange(ts_dens_EJ + 
                                        ggtitle("No transformation"), 
                                    ts_dens_bc_EJ + ggtitle("Box-Cox transformation"), 
                                    ncol = 2)

dense_cow <- cowplot::plot_grid(dense_all, 
                   dense_EJ, 
                   ncol = 1, 
                   labels = c('All', 'EJ'),
                   label_colour = "#004c93",
                   rel_heights = c(1, 1))

save(dense_cow, 
     file = here::here('09_simulation_and_approximation-cdf/01_paper/density_plots.Rdata'))



