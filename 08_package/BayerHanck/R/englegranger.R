englegranger <- function(formula, data, lags){

  eg_lm <- lm(formula, data = data)
  eg_res <- eg_lm$residuals
  eg_adf <- urca::ur.df(eg_res, lags = lags)
  return(as.numeric(eg_adf@teststat))

}
