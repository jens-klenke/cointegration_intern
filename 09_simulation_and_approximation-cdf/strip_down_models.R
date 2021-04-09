Data <- readRDS(here::here('09_simulation_and_approximation-cdf/Data.rds'))

# split dataset in cases
data_case_1 <- Data %>%
    dplyr::filter(case == 1)

data_case_2 <- Data %>%
    dplyr::filter(case == 2)

data_case_3 <- Data %>%
    dplyr::filter(case == 3)

### test try data set
test_data <- Data%>%
    dplyr::sample_n(10)

#### models ####
fit_case_1 <- lm(p_value_E_G ~ poly(stat_E_G, 3) + poly(k, 3) ,  data = data_case_1)
summary(fit_case_1)

### predict
predict(fit_case_1, test_data)

##### linear / poly / lm model 

clean_lm <- function(object) {
    object$y = c()
    object$model = c()
    
    object$residuals = c()
    object$fitted.values = c()
    object$effects = c()
    object$qr$qr = c()  
    object$linear.predictors = c()
    object$weights = c()
    object$prior.weights = c()
    object$data = c()
    
    
    object$family$variance = c()
    object$family$dev.resids = c()
    object$family$aic = c()
    object$family$validmu = c()
    object$family$simulate = c()
    attr(object$terms,".Environment") = c()
    attr(object$formula,".Environment") = c()
    
    object
}

#### GAM models ####
fit_gam_case_1 <- gam::gam(p_value_E_G ~ ns(stat_E_G, 3),  data = data_case_1)

fit_gam_case_1 <- gam_case.1_bc_poly_log_m_10_poly_m_sqrt

fit_gam_case_1$fitted.values <- NULL
fit_gam_case_1$residuals <- NULL
fit_gam_case_1$effects <- NULL
fit_gam_case_1$weights <- NULL
fit_gam_case_1$weights <- NULL
fit_gam_case_1$rank <- NULL
fit_gam_case_1$assign <- NULL
fit_gam_case_1$qr$qr <- NULL
fit_gam_case_1$df.residual <- NULL
fit_gam_case_1$additive.predictors <- NULL
fit_gam_case_1$y <- NULL
fit_gam_case_1$family <- NULL
fit_gam_case_1$prior.weights <- NULL
fit_gam_case_1$model <- NULL
fit_gam_case_1$data <- NULL

object_size(fit_gam_case_1)





