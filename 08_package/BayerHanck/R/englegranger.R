#' Engle-Granger Test
#'
#' Executes Engle-Granger Test.
#'
#' @param formula An object of class "formula" to describe the model.
#' @param data An optional data frame containing the variables in the model.
#' @param lags Number of lags to be included.
#' @param trend Type of deterministic component to be inlcuded, "none" for no deterministics,
#' "const" for a constant and "trend" for a constant plus trend.
#'
#' @return \code{englegranger} returns an object of class "co.test".
#' @export
#'
#' @references Engle, R. and Granger, C. (1987), Co-integration and Error Correction: Representation, Estimation, and Testing, Econometrica 55(2), 251-76.
#'
#' @examples englegranger(linvestment ~ lincome + lconsumption, data = Lutkepohl)
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

  out <- list(test.stat = test.stat,
              lags = lags,
              trend = trend,
              test = "Engle-Granger")
  class(out) <- c("co.test", "list")
  cat(c("----------------------------------------------------------",
        "Engle-Granger Test",
        "----------------------------------------------------------",
        paste(c("Value of test statistic:", round(test.stat, 4)), collapse = " ")),
        sep = "\n")
  invisible(out)
}
