#-- functions ----
# Design matrices
model.frame.fastLm <- function (formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 
                        0)]
    fcall <- formula$call
    m <- match(c("formula", "data", "subset", "weights", 
                 "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- quote(stats::model.frame)
    fcall$xlev <- formula$xlevels
    fcall$formula <- terms(formula)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env)) 
        env <- parent.frame()
    eval(fcall, env)
}

model.matrix.fastLm <- function (object, ...) {
    data <- model.frame.fastLm(object, ...)
    dots <- list(...)
    dots$data <- dots$contrasts.arg <- NULL
    do.call("model.matrix.default", c(list(object = object, 
                                           data = data, contrasts.arg = list(k_dummy = "contr.treatment")), 
                                      dots))
}

# inverse BoxCox function
invBoxCox <- function(x){
    x <- if (lambda_p == 0) exp(as.complex(x)) else (lambda_p*as.complex(x) + 1)^(1/lambda_p)
    return(Re(x))
}

# get lambda 
get_lambda <- function(data, case_w, art, test_w){
    data %>%
        dplyr::filter(case == case_w, 
                      test.type == test_w) %>%
        dplyr::select(all_of(art)) %>%
        dplyr::pull()
}


get_p_trans <- function(model){
    model %>%
        purrr::pluck('formula') %>%
        purrr::pluck(2)
}

get_model <- function(trendtype, test){
    models %>%
        dplyr::filter(case == trendtype, 
                      test.type == test) %>%
        dplyr::select(model) %>%
        dplyr::pull() %>%
        purrr::pluck(1)
}


get_p_value <- function(bh.test, trendtype, test.type, k){
    # trendtype = case
    
    # getting lambda for p_values
    lambda_p <- get_lambda(models, trendtype, 'p', 'all')
    # getting lambda for stat
    lambda_stat <- get_lambda(models, trendtype, 'stat', test.type)
    # saving model 
    model <- get_model(trendtype, test.type) 
    # dependent var 
    dep_var <- get_p_trans(model)
    
    # generating data set 
    new_data <- tibble(p_value_Fisher = 1L, 
                       stat_Fisher_all = bh.test, 
                       stat_Fisher_E_J = bh.test, 
                       k = k) %>%
        dplyr::mutate(k_dummy = as.factor(k),
                      stat_Fisher_E_J_bc = ((stat_Fisher_E_J^lambda_stat)-1)/lambda_stat,
                      stat_Fisher_all_bc = ((stat_Fisher_all^lambda_stat)-1)/lambda_stat)
    # approximation of the model 
    p.value_raw <- as.vector(model.matrix.fastLm(object = model, data = new_data) %*% coef(model))
    
    p.value_trans <- if(stringr::str_detect(dep_var, '_lg')) {
        exp(p.value_raw)
    } else {
        if(stringr::str_detect(dep_var, '_bc')){
            invBoxCox(p.value_raw)
        } else {p.value_raw}
    }
    
    p.value <- ifelse(p.value_trans >= 1, 9.9999e-1, ifelse(p.value_trans <= 0, 1e-12, p.value_trans))
    return(p.value)
}

bh.test <- 10 
trendtype <- 1
test.type <- 'all'
k <- 5

get_p_value(10, 1, 'all', 5)





