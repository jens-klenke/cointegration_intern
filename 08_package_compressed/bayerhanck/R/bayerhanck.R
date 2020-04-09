#' Joint Test-Statistic for the Null of Non-Cointegration
#'
#' Produces a joint test-statistic for the null of non-cointegration, aggregating
#' various cointegration tests.
#'
#' @param formula An object of class \code{\link[stats]{formula}} to describe the model.
#' @param data An optional data frame containing the variables in the model.
#' @param lags Number of lags to be included.
#' @param trend Type of deterministic component to be included, "none" for no deterministic,
#' "const" for a constant and "trend" for a constant plus trend.
#' @param test Selection of tests to choose from. Choices are either "eg-j", for \code{\link{englegranger}}
#' and \code{\link{johansen}}, or "all", for \code{\link{englegranger}}, \code{\link{johansen}},
#' \code{\link{banerjee}} and \code{\link{boswijk}}.
#' @param crit Level for the critical value of the test to be reported.
#'
#' @return \code{bayerhanck} returns an object of classes \code{"bh.test"} and \code{"list"}.
#'
#' The function \code{summary} is used to print a summary, whereas the cumulative distribution
#' under the null hypothesis can be plotted with \code{plot}.
#'
#' @export
#'
#' @references Bayer, C. and Hanck, C. (2009), Combining Non-Cointegration tests, METEOR RM 09/012, University of Maastricht.
#'
#' @examples
#' data("mts-examples", package="MTS")
#' bayerhanck(sp ~ ibm + ko, data = ibmspko)
bayerhanck <- function(formula, data, lags = 1, trend = "const", test = "all", crit = 0.05) {

  #-----------------------------------------------------------------------------------------
  # Check Syntax
  #-----------------------------------------------------------------------------------------
  mf <- match.call()
  m <- match(c("formula", "data"), names(mf), 0L)
  if (is.null(data))
    stop()
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  mf <- na.omit(mf)
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)[, -1]
  nvar <- ncol(cbind(y, x))
  trend <- match.arg(trend,
                     choices = c("none", "const", "trend"))
  test <- match.arg(test,
                    choices = c("eg-j", "all"))
  crit <- match.arg(as.character(crit),
                    c(0.01, 0.05, 0.10))

  lag <- lags
  if (lag < 0)
    stop("Lags must be set to a non negative value.")
  if (crit < 0)
    stop("Level of critical value must be set to a non negative.")

    #-----------------------------------------------------------------------------------------
  # Code trendtypes
  #-----------------------------------------------------------------------------------------
  if (identical(trend, "none")){
    ending = -1
    trendtype = 1
  } else if (identical(trend, "const")){
    ending = ""
    trendtype = 2
  } else if (identical(trend, "trend")){
    ending = "FORMULAR"
    trendtype = 3
  }

  #-----------------------------------------------------------------------------------------
  # Call Tests
  #-----------------------------------------------------------------------------------------
  test.stat <- rep(NA, 4)
  names(test.stat) <- c("Engle-Granger", "Johansen", "Banerjee", "Boswijk")

  invisible(capture.output(
  if (identical(test, "eg-j"))
    test.stat[1:2] <- c(englegranger(formula = formula, data = data, lags = lags, trend = trend)$test.stat,
                        johansen(formula = formula, data = data, lags = lags, trend = trend)$test.stat)
  ))
  invisible(capture.output(
  if (identical(test, "all"))
    test.stat[1:4] <- c(englegranger(formula = formula, data = data, lags = lags, trend = trend)$test.stat,
                        johansen(formula = formula, data = data, lags = lags, trend = trend)$test.stat,
                        banerjee(formula = formula, data = data, lags = lags, trend = trend)$test.stat,
                        boswijk(formula = formula, data = data, lags = lags, trend = trend)$test.stat)
  ))

  test.stat <- test.stat
  pval.stat <- test.stat

  #-----------------------------------------------------------------------------------------
  # Obtain P-Values
  #-----------------------------------------------------------------------------------------

  N <- nrow(null_dist)

  basecase <- 44 * (trendtype - 1) + 4 * (nvar - 2)

  for (i in 1:4) {
    case = basecase + i
    if (i %in% c(1, 3)) {
      n <- sum(pval.stat[i] > null_dist[, case])
      pval.stat[i] <-  (n/N) + .000000000001
    } else {
      if (i %in% c(2, 4)) {
        n <- sum(pval.stat[i] < null_dist[, case])
        pval.stat[i] <-  (n/N) + .000000000001
      }
    }
  }

  #-----------------------------------------------------------------------------------------
  # Calculate Bayer-Hanck Fisher Statistics
  #-----------------------------------------------------------------------------------------

  #### Select matrices of critical value ####

  if (identical(crit, 0.01)){
    crit_val_1 <- crit_val_1_0.01
    crit_val_2 <- crit_val_2_0.01
  }

  if (identical(crit, 0.05)){
    crit_val_1 <- crit_val_1_0.05
    crit_val_2 <- crit_val_2_0.05
  }

  if (identical(crit, 0.10)){
    crit_val_1 <- crit_val_1_0.10
    crit_val_2 <- crit_val_2_0.10
  }

  #### compute statistics ####
  if (identical(test, "eg-j"))
    bh.test <- -2*sum(log(pval.stat[1:2]))
  if (identical(test, "all"))
    bh.test <- -2*sum(log(pval.stat[1:4]))

  # degrees of freedom
  n_var <- nvar - 1

  ### obtain critical Value
  if (identical(test, "eg-j"))
    crit.val <- crit_val_1[n_var, trendtype]
  if (identical(test, "all"))
    crit.val <- crit_val_2[n_var, trendtype]

  #-----------------------------------------------------------------------------------------
  # Display Results
  #-----------------------------------------------------------------------------------------
  out <- list(bh.test = bh.test,
              test.stat = test.stat[complete.cases(test.stat)],
              pval.stat = pval.stat[complete.cases(test.stat)],
              crit.val = crit.val,
              formula = formula,
              lags = lags,
              trend = trend,
              crit = crit,
              nvar = nvar,
              test.type = test,
              K = n_var,
              bh.crit.val = crit.val,
              basecase = basecase)
  class(out) <- c("bh.test", "list")
  cat(c("----------------------------------------------------------",
        "Bayer-Hanck Test for Non-Cointegration",
        "----------------------------------------------------------",
        paste(c("Value of the Fisher Type Test statistic:", round(bh.test, 4)),
              collapse = " "),
      paste(c("Critical" ,round(crit*100, 0) ,"% Value of the Fisher Type Test:", round(crit.val, 4)),
            collapse = " ")),
      sep = "\n")
  invisible(out)
}
