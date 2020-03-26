bayerhanck <- function(formula, data, trend, lags, test, crit) {
  mf <- model.frame(formula)
  y <- model.response(mf)

  # Engle Granger Test

  # Johansen Test

  # Boswijk/Banerjee Test

}


y <- 1:4
x <- cbind(5:8, 10:13)
z <- cbind(1:4, 20:23)

form <- y ~ x + z
mf <- model.frame(form)
y_mf <- model.response(mf)
?model.response
