bayerhanck <- function(formula, data, trend, lags, test, crit, ...) {

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
  test.stat <- as.numeric()

  if (test == "banerjee")
    test.stat <- banerjee(formula = formula, data = data, lags = lags, trend = trend)
  if (test == "boswijk")
    test.stat <- boswijk(formula = formula, data = data, lags = lags, trend = trend)
  if (test == "englegranger")
    test.stat <- englegranger(formula = formula, data = data, lags = lags, trend = trend)
  if (test == "johansen")
    test.stat <- johansen(x = model.frame(formula = formula, data = data),
                          lags = lags, trend = trend)

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

bayerhanck(data = df, formula = linvestment ~ lincome + lconsumption,
           lags = 1, test = "johansen")

