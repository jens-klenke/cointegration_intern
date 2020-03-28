bayerhanck <- function(formula, data, trend = "const", lags = 1, test = "all", crit) {

  #-----------------------------------------------------------------------------------------
  # Check Syntax
  #-----------------------------------------------------------------------------------------
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

  #-----------------------------------------------------------------------------------------
  # Call Tests
  #-----------------------------------------------------------------------------------------
  test.stat <- rep(NA, 4)
  names(test.stat) <- c("englegranger", "johansen", "banerjee", "boswijk")

  if ("englegranger" %in% test)
    test.stat[1] <- englegranger(formula = formula, data = data, lags = lags, trend = trend)[[1]]
  if ("johansen" %in% test)
    test.stat[2] <- johansen(formula = formula, data = data, lags = lags, trend = trend)[[1]]
  if ("banerjee" %in% test)
    test.stat[3] <- banerjee(formula = formula, data = data, lags = lags, trend = trend)[[1]]
  if ("boswijk" %in% test)
    test.stat[4] <- boswijk(formula = formula, data = data, lags = lags, trend = trend)[[1]]
  if (identical(test, "all"))
    test.stat[1:4] <- c(englegranger(formula = formula, data = data, lags = lags, trend = trend)[[1]],
                        johansen(formula = formula, data = data, lags = lags, trend = trend)[[1]],
                        banerjee(formula = formula, data = data, lags = lags, trend = trend)[[1]],
                        boswijk(formula = formula, data = data, lags = lags, trend = trend))[[1]]

  #-----------------------------------------------------------------------------------------
  # Obtain P-Values
  #-----------------------------------------------------------------------------------------
  load(here::here("/null_dist.rda"))
  crit_val <- rmatio::read.mat(here::here("/critical_values.mat"))

  pval <- test.stat
  basecase <- 44 * (trendtype - 1) + 4 * (nvar - 2)

  for (i in 1:4) {
    case = basecase + i
    if (i %in% c(1, 3)) {
      n <- sum(stat[i] > null_dist[, case])
      pval[i] <-  (n/N) + .000000000001
    } else {
      if (i %in% c(2, 4)) {
        n <- sum(stat[i] < null_dist[, case])
        pval[i] <-  (n/N) + .000000000001
      }
    }
  }

  #-----------------------------------------------------------------------------------------
  # Calculate Bayer-Hanck Fisher Statistics
  #-----------------------------------------------------------------------------------------

  #-----------------------------------------------------------------------------------------
  # Display Results
  #-----------------------------------------------------------------------------------------
  list(test.stat = test.stat,
       p.val = pval)

}


bayerhanck(linvestment ~ lincome + lconsumption, data = df, test = c("englegranger",
                                                                           "banerjee"))
