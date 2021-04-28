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
k = 2 
bh.test <- a$bh.test
trendtype  <- 1
test.type <- a$test.type


get_p_value <- function(bh.test, trendtype, test.type){
    #lambda_p <<- get_lambda(trendtype, 'p', 'all')
    
    new_data <- tibble(stat_Fisher_all = bh.test, 
                   stat_Fisher_E_J = bh.test, 
                   k =  k ) %>%
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

#---- Test ----
get_p_value(20.1, 1, 'all')
get_p_value(50, 1, 'all')
get_p_value(55, 1, 'all')
get_p_value(100, 1, 'all')
get_p_value(200, 1, 'all')

Test_stat <- 1:400
p_values_1_all <- tibble(
    case = 1,
    p_value = purrr::map_dbl(Test_stat, ~get_p_value(., 1, 'all')), 
    test_stat = Test_stat)
p_values_2_all <- tibble(
    case = 2,
    p_value = purrr::map_dbl(Test_stat, ~get_p_value(., 2, 'all')), 
    test_stat = Test_stat)
p_values_3_all <- tibble(
    case = 3,
    p_value = purrr::map_dbl(Test_stat, ~get_p_value(., 3, 'all')), 
    test_stat = Test_stat)

p_values_all <- rbind(
    p_values_1_all, p_values_2_all, p_values_3_all
)

p_values_1_e_j <- tibble(
    case = 1,
    p_value = purrr::map_dbl(Test_stat, ~get_p_value(., 1, 'e_j')), 
    test_stat = Test_stat)
p_values_2_e_j <- tibble(
    case = 2,
    p_value = purrr::map_dbl(Test_stat, ~get_p_value(., 2, 'e_j')), 
    test_stat = Test_stat)
p_values_3_e_j <- tibble(
    case = 3,
    p_value = purrr::map_dbl(Test_stat, ~get_p_value(., 3, 'e_j')), 
    test_stat = Test_stat)

p_values_e_j <- rbind(
    p_values_1_e_j, p_values_2_e_j, p_values_3_e_j
)

test_plot <- function(p_values){
    p_values %>%
        ggplot(aes(x = test_stat, y = p_value)) +
        geom_line(color = '#004c93') +
        labs(x = 'Test Statistic', y = 'Approximated p-values \n')+
        theme_bw()+
        facet_wrap(~case)+
        theme(panel.spacing = unit(1, "lines"),
              strip.background = element_rect(colour = 'black',
                                              fill = '#004c93'),
              strip.text.x = element_text(size = 12, color = 'white' # , face = "bold.italic"
              ))
}

test_plot(p_values_all)
test_plot(p_values_e_j)



## Janine ausführen und speichern 

# Load Simulation Data 
if(Sys.info()['nodename'] == "DELL-ARBEIT") { # Jens 
    Data <- readRDS('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_100k.rds')
    # load('C:\\Users\\Jens-\\Dropbox\\jens\\BayerHanck\\Data_1_m.RData')
} else if(Sys.info()['user'] == "Janine") { # Janine
    load("/Users/Janine/Desktop/BayerHanck/Data_1_m.RData")
} else if(Sys.info()['nodename'] == "OEK-TS01") { # Server
    load('D:\\Klenke\\Data_1_m.RData')
}

crit_val_all <- Data %>%
    dplyr::group_by(case, k) %>%
    dplyr::select(stat_Fisher_all) %>%
    dplyr::slice_max(stat_Fisher_all , n = 100) %>%
    dplyr::filter(stat_Fisher_all == min(stat_Fisher_all))

crit_val_e_j <- Data %>%
    dplyr::group_by(case, k) %>%
    dplyr::select(stat_Fisher_E_J) %>%
    dplyr::slice_max(stat_Fisher_E_J , n = 100) %>%
    dplyr::filter(stat_Fisher_E_J == min(stat_Fisher_E_J))

