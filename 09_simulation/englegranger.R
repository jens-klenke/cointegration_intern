englegranger <- function(formula, data, lags = 1, trend = "const"){
  eg_lm <- lm(formula, data = data)
  eg_res <- eg_lm$residuals
  eg_adf <- urca::ur.df(eg_res, lags = lags)
  test.stat <- as.numeric(eg_adf@teststat)
  names(test.stat) <- "englegranger"

  list(test.stat = test.stat)
}

