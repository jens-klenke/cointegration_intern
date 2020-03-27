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
  test.stat <- numeric()
  if (test == "all")
    test.stat[1:4] <- c(englegranger(formula = formula, data = data, lags = lags, trend = trend),
                        #johansen(formula = formula, data = data, lags = lags, trend = trend),
                        banerjee(formula = formula, data = data, lags = lags, trend = trend),
                        boswijk(formula = formula, data = data, lags = lags, trend = trend))
  if (test == "banerjee")
    test.stat[1] <- banerjee(formula = formula, data = data, lags = lags, trend = trend)
  if (test == "boswijk")
    test.stat[1] <- boswijk(formula = formula, data = data, lags = lags, trend = trend)
  if (test == "englegranger")
    test.stat[1] <- englegranger(formula = formula, data = data, lags = lags, trend = trend)
  #if (test == "johansen")
  #  test.stat[1] <- johansen(formula = formula, data = data, lags = lags, trend = trend)

  #-----------------------------------------------------------------------------------------
  # Obtain P-Values
  #-----------------------------------------------------------------------------------------

  #-----------------------------------------------------------------------------------------
  # Calculate Bayer-Hanck Fisher Statistics
  #-----------------------------------------------------------------------------------------

  #-----------------------------------------------------------------------------------------
  # Display Results
  #-----------------------------------------------------------------------------------------
}

bayerhanck(linvestment ~ lincome + lconsumption, data = df, test = "all")

