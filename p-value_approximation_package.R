#---- Preliminary ---- 

# packages
source(here::here('01_code/packages/packages.R'))

# load models 
load(here::here("09_simulation_and_approximation-cdf/p_value_approx.RData"))

# new data 

# bayerhanck names 
## trendtype  == case
## stat = bh.test
## k = nvar
## test "eg-j" , "all"

trendtype <- 2
case <- 2
test = 'eg-j'
nvar = 1 
stat <- 0.5

#-- simulate own data from package output ----
new_data <- tibble(
    k = nvar,
    stat_Fisher_all = NA,
    stat_Fisher_all_lg = NA,
    stat_Fisher_all_bc = NA, 
    stat_Fisher_E_J = NA, 
    stat_Fisher_E_J_lg = NA,
    stat_Fisher_E_J_bc = NA)%>%
    dplyr::mutate(stat_Fisher_E_J = ifelse(identical(test, 'eg-j'), stat, NA),
                  stat_Fisher_all = ifelse(identical(test, 'all'), stat, NA)
                  )

# bc und log berechnen!! 
if(identical(test, 'eg-j') &  identical(trendtype, 1)){
    p.value <- predict(mod_E_J_case_1, new_data)
}


function(test, case, )


if(identical(test, 'eg-j') &  identical(trendtype, 1)){
    p.value <- predict(mod_E_J_case_1, new_data)
}



    




