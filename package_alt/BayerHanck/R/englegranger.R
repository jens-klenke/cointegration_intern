englegranger <- function(formula, data, lags = 1, trend = "const"){

  #-----------------------------------------------------------------------------------------
  # Trend Specification
  #-----------------------------------------------------------------------------------------
  if (identical(trend, "none")) {
    eg_lm <- lm(update(formula, ~. -1), data = data, na.action = na.omit)
  } else if (identical(trend, "const")) {
    eg_lm <- lm(formula, data = data, na.action = na.omit)
  } #else if (identical(trend, "trend")) {}

  #-----------------------------------------------------------------------------------------
  # Engle Granger Test
  #-----------------------------------------------------------------------------------------
  eg_res <- eg_lm$residuals
  eg_adf <- urca::ur.df(eg_res, lags = lags)
  test.stat <- as.numeric(eg_adf@teststat)
  names(test.stat) <- "englegranger"

  list(test.stat = test.stat)
}



