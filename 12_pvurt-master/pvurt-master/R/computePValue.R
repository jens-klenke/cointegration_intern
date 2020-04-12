#' P-value
#'
#' Function to approximate the p-value for augmented Dickey-Fuller test.
#'
#' Based on the chosen model (GAM or polynomial), the function returns the
#' approximated p-value. Default is GAM model.
#'
#' @name computePValue
#'
#' @rdname computePValue
#'
#' @param object Numeric value or an object (\code{fHTEST}, \code{ur.df} or
#'   \code{htest}) for which p-value needs to approximated.
#' @param ... Further arguments passed to methods.
#' @param n Sample size.
#' @param type The type of unit root test. Currently supports: \code{nc} for
#'   test without drift and trend, \code{c} for test with only drift and
#'   \code{ct} for test with both drift and trend.
#' @param model The model type to be used for approximation. Available is GAM
#'   and polynomial regression. If \code{gam} is chosen, then \code{d} has no
#'   effect.
#' @param d The degree for polynomial. \code{d} must be \eqn{\ge 3} and \eqn{\le
#'   6}. If \code{gam} is chosen, then \code{d} has no effect.
#'
#' @examples
#' library(pvurt)
#' y <- arima.sim(model = list(order = c(0, 1, 0)), n = 100)
#'
#' # Test type: with drift and trend
#' # package: fUnitRoots
#' library(fUnitRoots)
#' computePValue(adfTest(y, lags = 3, type = "ct"))
#' computePValue(unitrootTest(y, lags = 3, type = "ct"))
#'
#' # package: urca
#' library(urca)
#' computePValue(ur.df(y, lags = 3, type = "trend"))
#' # print summary
#' summary(computePValue(ur.df(y, lags = 3, type = "trend")))
#'
#' # package: tseries
#' library(tseries)
#' computePValue(adf.test(y, alternative = "stationary", k = 3))
#'
#' # no packages
#' tStat <- -2.239
#' sampleSize <- 100
#' computePValue(tStat, n = sampleSize, model = "gam", type = "ct")


NULL

#' @rdname computePValue
#' @export
computePValue <- function(object, ...){
  UseMethod("computePValue")
}

#' @rdname computePValue
#' @method computePValue numeric
#' @export
computePValue.numeric <- function(
  object, n, type = c("nc", "c", "ct"), model = c("gam", "poly"), d = NULL, ...
  ){
  type <- match.arg(type)
  model <- match.arg(model)
  d <- d
  p_ <- computePValueMain(object = object, n = n, model = model, type = type, d = d)[[1]]
  return(p_)
}

#' @rdname computePValue
#' @importFrom stringr str_match
#' @import timeSeries
#' @export
computePValue.fHTEST <- function(
  object, model = c("gam", "poly"), type = c("nc", "c", "ct"), d = NULL, ...
){
  attrObject <- attributes(object)
  model <- match.arg(model)

  x <- attrObject$data$x
  y <- diff(x)
  n <- length(y)

  typeCall <- stringr::str_match(deparse(attrObject$call), 'type = \\"([nct]*)\\"')[1,2]

  if(!(is.na(typeCall))){
    if(typeCall == "ctt") {
      type_ <- "ct"
      warning("'ctt' not supported. Using 'ct' instead.")
    } else {
      type_ <- typeCall
    }
  } else {
    type_ <- match.arg(type)
    warning(paste(
      "Test type could not be determined. Specify using the 'type' argument.",
      "Currently supported: nc, c, ct. Uses 'nc' (Without drift and trend) by default."
    ))
  }

  t <- attrObject$test$statistic[[1]]

  p_ <- computePValueMain(object = t, n = n, model = model, type = type_, d = d)[[1]]

  oldNames <- names(attrObject$test$p.value)
  if(oldNames[1] == "") oldNames[1] = "t"
  newNames <- c(oldNames, "pvurt")

  oldPVAL <- attrObject$test$p.value
  newPVAL <- c(oldPVAL, p_)
  names(newPVAL) <- newNames

  attrObject$test$p.value <- newPVAL

  new("fHTEST", data = list(x = attrObject$data$x), test = attrObject$test,
      title = as.character(attrObject$title), description = description())
}

#' @rdname computePValue
#' @export
computePValue.ur.df <- function(
  object, model = c("gam", "poly"), d = NULL, ...
){

  attrObject <- attributes(object)
  model <- match.arg(model)

  y <- attrObject$y
  z <- diff(y)
  n <- length(z)

  typeVector <- c("nc", "c", "ct")
  names(typeVector) <- c("none", "drift", "trend")

  type <- typeVector[[attrObject$model]]
  t <- attrObject$teststat[[1]]

  p_ <- computePValueMain(object=t, n = n, model = model, type = type, d = d)[[1]]

  pValueVector <- rep(NA, ncol(attrObject$teststat))
  pValueVector[1] <- p_
  pValueVector <- as.matrix(pValueVector)

  new("ur.df.pvurt", y = attrObject$y, model = attrObject$model, cval = attrObject$cval,
      lags = attrObject$lags, teststat = attrObject$teststat, approxpval = pValueVector,
      testreg = attrObject$testreg, res = attrObject$res,
      test.name = "Augmented Dickey-Fuller Test")
}

#' @rdname computePValue
#' @export
computePValue.htest <- function(
  object, model = c("gam", "poly"), d = NULL, ...
){
  t <- object["statistic"][[1]][[1]]
  type <- "ct"
  model <- match.arg(model)

  x <- eval(parse(text = object["data.name"][[1]][[1]]))
  n <- length(diff(x))

  p_ <- computePValueMain(object = t, n = n, model = model, type = type, d = d)[[1]]

  alternative <- object["alternative"][[1]][[1]]
  if(alternative == "stationary"){
    pvurtValue <- p_
  } else if(alternative == "explosive"){
    pvurtValue <- 1 - p_
  }

  STAT <- object["statistic"][[1]][[1]]
  PARAMETER <- object["parameter"][[1]][[1]]
  PVAL <- object["p.value"][[1]][[1]]
  METHOD <- object["method"][[1]][[1]]
  DNAME <- object["data.name"][[1]][[1]]
  names(STAT) <- "Dickey-Fuller"
  names(PARAMETER) <- "Lag order"

  new("htest.pvurt", statistic = STAT, parameter = PARAMETER, alternative = alternative,
      p.value = PVAL, pvurt = pvurtValue, method = METHOD, data.name = DNAME)
}


