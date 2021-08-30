# packages 
source(here::here('01_code/packages/packages.R'))

# loading data 
load(file = here::here('09_simulation_and_approximation-cdf/plot_data_correctur.RData'))

# load model matrix
load(here::here('09_simulation_and_approximation-cdf/models.RData'))

# load sub-functions
source(here::here('09_simulation_and_approximation-cdf/functions_package.R'))

# load graphic functions
source(here::here('09_simulation_and_approximation-cdf/functions_graphic.R'))

# parallel 
plan(multisession, workers = 4)

#-- plot sim-p-values vs predicted p-values ----

# diagnose plots
p.sim_p.aprox_all <- own_plot(data_all)
p.sim_p.aprox_all_0.2 <- own_plot_0.2(data_all)

p.sim_p.aprox_all.k1 <- own_plot_k1(data_all)
p.sim_p.aprox_all_0.2.k1 <- own_plot_0.2_k1(data_all)

# diagnose plots
p.sim_p.aprox_e_j <- own_plot(data_e_j)
p.sim_p.aprox_e_j_0.2 <- own_plot_0.2(data_e_j)

p.sim_p.aprox_e_j.k1 <- own_plot_k1(data_e_j)
p.sim_p.aprox_e_j_0.2.k1 <- own_plot_0.2_k1(data_e_j)

#-- plot stat vs predicted p-values ----

# generating data sets
p_value_plots_e_j <- p_value_test_fun('E_J')
p_value_plots_all <- p_value_test_fun('all')

# plot
plot_p_stat_e_j <- plot_p_stat(p_value_plots_e_j)
plot_p_stat_all <- plot_p_stat(p_value_plots_all)

plot_p_stat_e_j_k.2 <- plot_p_stat_k.2(p_value_plots_e_j)
plot_p_stat_all_k.2 <- plot_p_stat_k.2(p_value_plots_all)

#-- saving plots ----
save(p.sim_p.aprox_all.k1, p.sim_p.aprox_all_0.2.k1, p.sim_p.aprox_e_j.k1, p.sim_p.aprox_e_j_0.2.k1, 
     p.sim_p.aprox_all, p.sim_p.aprox_all_0.2, p.sim_p.aprox_e_j, p.sim_p.aprox_e_j_0.2,
     plot_p_stat_e_j, plot_p_stat_all, plot_p_stat_e_j_k.2, plot_p_stat_all_k.2,
     file = here::here('09_simulation_and_approximation-cdf/01_paper/paper_plots.RData'))
