# packages 
source(here::here('01_code/packages/packages.R'))

# loading data 
load(file = here::here('09_simulation_and_approximation-cdf/01_paper/paper_plots.RData'))

#-- p-value: sim vs approx ----

# p.sim-approx all
tikz(here::here('09_simulation_and_approximation-cdf/01_paper/includes/p.sim_p.aprox_all.tex'), width = 10, height = 7.5)
print(p.sim_p.aprox_all)
dev.off()

# p.sim-approx all 0.2
tikz(here::here('09_simulation_and_approximation-cdf/01_paper/includes/p.sim_p.aprox_all_0.2.tex'), width = 10, height = 10)
print(p.sim_p.aprox_all_0.2)
dev.off()
