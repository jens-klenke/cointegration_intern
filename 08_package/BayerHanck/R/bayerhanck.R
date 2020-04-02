#' Joint Test-Statistic for the Null of Non-Cointegration
#'
#' Produces a joint test-statistic for the null of non-cointegration, aggregating
#' various cointegration tests.
#'
#' @param formula An object of class "formula" to describe the model.
#' @param data An optional data frame containing the variables in the model.
#' @param lags Number of lags to be included.
#' @param trend Type of deterministic component to be inlcuded, "none" for no deterministics,
#' "const" for a constant and "trend" for a constant plus trend.
#' @param test Selection of tests to choose from.
#' @param crit Level for the critical value of the test to be reported.
#'
#' @return
#' @export
#'
#' @references Bayer, C. and Hanck, C. (2009), Combining Non-Cointegration tests, METEOR RM 09/012, University of Maastricht.
#'
#' @examples
bayerhanck <- function(formula, data, lags = 1, trend = "const", test = "all", crit = 0.05) {

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
  x <- model.matrix(mt, mf)[, -1]
  nvar <- ncol(cbind(y, x))

  #if (nrow(x) == 0L)
  #  stop("0 (non-NA) cases")
  #if (NROW(y) != nrow(x))
  #  stop("Incompatible dimensions")
  #if (trend != "none" & trend != "constant" & trend != "trend")
  #  stop("Trend cannot be specified")
  #if (lags < 1)
  #  stop("Lags must be >= 1")
  #if (lags < nrow(x))
  #  stop("Lags cannot exceed number of observations")
  #if (crit < 0)
  #  stop("Level of critical value cannot be negative")

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
  names(test.stat) <- c("englegranger", "johansen", "banerjee", "boswijk")

  if ("englegranger" %in% test)
    test.stat[1] <- englegranger(formula = formula, data = data, lags = lags, trend = trend[[1]])
  if ("johansen" %in% test)
    test.stat[2] <- johansen(formula = formula, data = data, lags = lags, trend = trend[[1]])
  if ("banerjee" %in% test)
    test.stat[3] <- banerjee(formula = formula, data = data, lags = lags, trend = trend[[1]])
  if ("boswijk" %in% test)
    test.stat[4] <- boswijk(formula = formula, data = data, lags = lags, trend = trend[[1]])
  if (identical(test, "all"))
    test.stat[1:4] <- c(englegranger(formula = formula, data = data, lags = lags, trend = trend)[[1]],
                        johansen(formula = formula, data = data, lags = lags, trend = trend)[[1]],
                        banerjee(formula = formula, data = data, lags = lags, trend = trend)[[1]],
                        boswijk(formula = formula, data = data, lags = lags, trend = trend)[[1]])
  pval.stat <- test.stat[complete.cases(test.stat)]
  test.stat <- test.stat[complete.cases(test.stat)]
  print(test.stat)

  #-----------------------------------------------------------------------------------------
  # Obtain P-Values
  #-----------------------------------------------------------------------------------------
  load('null_dist.rda')

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

  #### Load critical values ####
  load('crit_values.rda')

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
  if (identical(test, "all"))
  b_h_stat_1 <- -2*sum(log(pval.stat[1:2]))
  b_h_stat_2 <- -2*sum(log(pval.stat[1:4]))

  # degrees of freedom
  nvar <- nvar - 1

  ### obtain critical Value
  cv_1 <- crit_val_1[nvar, trendtype]
  cv_2 <- crit_val_2[nvar, trendtype]

  #-----------------------------------------------------------------------------------------
  # Display Results
  #-----------------------------------------------------------------------------------------
  list(test.stat = test.stat,
       pval.stat = pval.stat)
  print(pval.stat)
}



