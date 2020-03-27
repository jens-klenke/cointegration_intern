johansen <- function(formula, data, type = "eigen", lags = 1, trend = "const"){

  #-----------------------------------------------------------------------------------------
  # Check Syntax
  #-----------------------------------------------------------------------------------------
  mf <- match.call()
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  offset <- model.offset(mf)
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)[, -1]
  x <- cbind(y, x)

  #-----------------------------------------------------------------------------------------
  # Johansen Test
  #-----------------------------------------------------------------------------------------
  jo_vec <- tsDyn::VECM(x, lag = lags, r = 2,
                        include = trend, # Bezeichnung Trend anpassen
                        estim = "ML")
  test.stat <- summary(tsDyn::rank.test(jo_vec, type = "eigen"))$eigen[1]
  print(test.stat)
}

