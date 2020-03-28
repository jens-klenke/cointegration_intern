#' Execute Engle-Granger Test
#'
#' @param formula An object of class "formula" to describe the model
#' @param data An optional data frame containing the variables in the model
#' @param lags Number of lags to be included
#' @param trend Type of deterministic component to be inlcuded, "none" for no deterministics,
#' "const" for a constant and "trend" for a constant plus trend
#'
#' @return
#' @export
#'
#' @examples
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



