#' Cointegration Tests
#'
#' Executes Cointegration Tests, which serve as underlying Tests for the Bayer-Hanck Fisher Type Statistic.
#'
#' @param formula An object of class \code{\link[stats]{formula}} to describe the model.
#' @param data An optional data frame containing the variables in the model.
#' @param lags Number of lags to be included.
#' @param trend Type of deterministic component to be inlcuded, "none" for no deterministics,
#' "const" for a constant and "trend" for a constant plus trend.
#' @param type Test to be conducted, either "eigen" or "trace".
#'
#' @return Returns an object of classes \code{"co.test"} and \code{"list"}.
#' @export
#'
#' @references Engle, R. and Granger, C. (1987), Co-integration and Error Correction: Representation, Estimation, and Testing, Econometrica 55(2), 251-76.
#'
#' Johansen, S. (1988), Statistical analysis of cointegration vectors, Journal of Economic Dynamics and Control 12(2-3), 231-254.
#'
#' Banerjee, A., Dolado, J. J. and Mestre, R. (1998), Error-correction Mechanism Tests for Cointegration in a Single-equation Framework, Journal of Times Series Analysis 19(3), 267-283.
#'
#' Boswijk, H. P. (1994), Testing for an unstable root in conditional and structural error correction models, Journal of Econometrics 63(1), 37-60.
#'
#' @examples
#' data("mts-examples", package="MTS")
#' englegranger(sp ~ ibm + ko, data = ibmspko)
englegranger <- function(formula, data, lags = 1, trend = "const"){

  #-----------------------------------------------------------------------------------------
  # Check Syntax
  #-----------------------------------------------------------------------------------------
  mf <- match.call()
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  mf <- na.omit(mf)
  y <- model.response(mf, "numeric")
  trend <- match.arg(trend,
                     choices = c("none", "const", "trend"))

  #-----------------------------------------------------------------------------------------
  # Trend Specification
  #-----------------------------------------------------------------------------------------
  if (identical(trend, "none")){
    eg_lm <- lm(update(formula, ~. -1), data = data, na.action = na.omit)
  } else if (identical(trend, "const")) {
    eg_lm <- lm(formula, data = data, na.action = na.omit)
  } else if (identical(trend, "trend")){
    tr <- seq_along(y)
    data$tr <- tr
    eg_lm <- lm(update(formula, ~. + tr), data = data, na.action = na.omit)
    data <- data[, -which(colnames(data) == "tr")]
  }

  #-----------------------------------------------------------------------------------------
  # Engle Granger Test
  #-----------------------------------------------------------------------------------------
  eg_res <- eg_lm$residuals
  eg_adf <- urca::ur.df(eg_res, lags = lags)
  test.stat <- as.numeric(eg_adf@teststat)

  out <- list(test.stat = test.stat,
              lags = lags,
              trend = trend,
              test = "Engle-Granger",
              formula = formula)
  class(out) <- c("co.test", "list")
  cat(c("----------------------------------------------------------",
        "Engle-Granger Test",
        "----------------------------------------------------------",
        paste(c("Value of test statistic:", round(test.stat, 4)), collapse = " ")),
        sep = "\n")
  invisible(out)
}
