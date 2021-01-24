#-- functions ----
bc_value <- function(test, case, variable){
    final_model_metrics %>%
        dplyr::filter(test == {{test}},
                      case == {{case}}) %>%
        dplyr::select({{variable}}) %>%
        pull()
    
}

# inverse BoxCox function
invBoxCox <- function(x, lambda){
  x <- if (lambda == 0) exp(as.complex(x)) else (lambda*as.complex(x) + 1)^(1/lambda)
  return(Re(x))
}

#-- Preliminary ---- 

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
    stat_Fisher_E_J_bc = NA) %>%
    dplyr::mutate(stat_Fisher_E_J = ifelse(identical(test, 'eg-j'), stat, NA),
                  stat_Fisher_all = ifelse(identical(test, 'all'), stat, NA)
                  )

#-- bc computing ----  

# eg-j test and case = 1
if(identical(test, 'eg-j') & identical(trendtype, 1)){
    new_data <- new_data %>%
      dplyr::mutate(stat_Fisher_E_J_bc = (((stat_Fisher_E_J^bc_value('eg-j', 1, 'bc_stat'))-1)/bc_value('eg-j', 1, 'bc_stat')),
                    stat_Fisher_E_J_lg = log(stat_Fisher_E_J)
                    )
}

# eg-j test and case = 2
if(identical(test, 'eg-j') &  identical(trendtype, 2)){
  new_data <- new_data%>%
    dplyr::mutate(stat_Fisher_E_J_bc = (((stat_Fisher_E_J^bc_value('eg-j', 2, 'bc_stat'))-1)/bc_value('eg-j', 2, 'bc_stat')),
                  stat_Fisher_E_J_lg = log(stat_Fisher_E_J)
    )
}

# eg-j test and case = 3
if(identical(test, 'eg-j') &  identical(trendtype, 3)){
  new_data <- new_data%>%
    dplyr::mutate(stat_Fisher_E_J_bc = (((stat_Fisher_E_J^bc_value('eg-j', 3, 'bc_stat'))-1)/bc_value('eg-j', 3, 'bc_stat')),
                  stat_Fisher_E_J_lg = log(stat_Fisher_E_J)
    )
}

# all test and case = 1
if(identical(test, 'all') & identical(trendtype, 1)){
  new_data <- new_data%>%
    dplyr::mutate(stat_Fisher_all_bc = (((stat_Fisher_all^bc_value('all', 1, 'bc_stat'))-1)/bc_value('all', 1, 'bc_stat')),
                  stat_Fisher_all_lg = log(stat_Fisher_all)
    )
}

# all test and case = 2
if(identical(test, 'all') & identical(trendtype, 2)){
  new_data <- new_data%>%
    dplyr::mutate(stat_Fisher_all_bc = (((stat_Fisher_all^bc_value('all', 2, 'bc_stat'))-1)/bc_value('all', 2, 'bc_stat')),
                  stat_Fisher_all_lg = log(stat_Fisher_all)
    )
}

# all test and case = 3
if(identical(test, 'all') & identical(trendtype, 3)){
  new_data <- new_data%>%
    dplyr::mutate(stat_Fisher_all_bc = (((stat_Fisher_all^bc_value('all', 3, 'bc_stat'))-1)/bc_value('all', 3, 'bc_stat')),
                  stat_Fisher_all_lg = log(stat_Fisher_all)
    )
}
    

#-- p - value ----

# eg-j and case = 1
if(identical(test, 'eg-j') &  identical(trendtype, 1)){
    p.value_lg <- predict(mod_E_J_case_1, new_data)
    
    p.value <- exp(p.value_lg)
}

# eg-j and case = 2
if(identical(test, 'eg-j') &  identical(trendtype, 2)){
  p.value_lg <- predict(mod_E_J_case_2, new_data)
  
  p.value <- exp(p.value_lg)
}

# eg-j and case = 3
if(identical(test, 'eg-j') &  identical(trendtype, 3)){
  p.value_lg <- predict(mod_E_J_case_3, new_data)
  
  p.value <- exp(p.value_lg)
}

# all and case = 1
if(identical(test, 'all') &  identical(trendtype, 1)){
  p.value_bc <- predict(mod_ALL_case_1, new_data)
  
  p.value <- invBoxCox(p.value_bc, bc_value('all', 3, 'bc_p'))
}

# all and case = 2
if(identical(test, 'all') &  identical(trendtype, 2)){
  p.value_lg <- predict(mod_ALL_case_2, new_data)
  
  p.value <- exp(p.value_lg)
}

# all and case = 3
if(identical(test, 'all') &  identical(trendtype, 3)){
  p.value_lg <- predict(mod_ALL_case_3, new_data)
  
  p.value <- exp(p.value_lg)
}


# p-value correction and 1- pvalue

p.value <- 1 - p.value

p.value <- ifelse(p.value <= 0, 1e-12, ifelse(p.value >= 1, 1-1e-06, p.value))