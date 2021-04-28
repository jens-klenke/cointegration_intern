# selecting lambda
get_lambda <- function(case_w, art, test_w){
    lambda_values %>%
        dplyr::filter(case == case_w, 
                      side == art, 
                      test == test_w) %>%
        dplyr::pull(value)
}

# getting model 
get_model <- function(case_w, art){
    models %>%
        dplyr::filter(case == case_w,
                      test == art) %>%
        dplyr::select(models) %>%
        dplyr::pull() %>%
        purrr::pluck(1)
}

# inverse BoxCox function
invBoxCox <- function(x, lambda = 0.7071139){
    x <- if (lambda == 0) exp(as.complex(x)) else (lambda*as.complex(x) + 1)^(1/lambda)
    return(Re(x))
}

# getting tansformation of p
get_p_trans <- function(case_w, art){
    models %>%
        dplyr::filter(test == art,
                      case == case_w) %>%
        dplyr::pull(response)
}
