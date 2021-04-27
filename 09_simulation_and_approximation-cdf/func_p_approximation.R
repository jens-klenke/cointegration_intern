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
invBoxCox <- function(x){
    x <- if (lambda_p == 0) exp(as.complex(x)) else (lambda_p*as.complex(x) + 1)^(1/lambda_p)
    return(Re(x))
}

# getting tansformation of p
get_p_trans <- function(case_w, art){
    models %>%
        dplyr::filter(test == art,
                      case == case_w) %>%
        dplyr::pull(response)
}
