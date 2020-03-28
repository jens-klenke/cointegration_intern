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
  jo_vec <- tsDyn::VECM(x, lag = lags,
                        include = trend, # Bezeichnung Trend anpassen
                        estim = "ML")
  jo_vec_sum <- summary(tsDyn::rank.test(jo_vec, type = type))
  test.stat <- as.numeric(jo_vec_sum$eigen[1])
  names(test.stat) <- "johansen"

  list(test.stat = test.stat,
       trace = cbind(jo_vec_sum$r, jo_vec_sum$trace, jo_vec_sum$trace_pval_T),
       eigen = cbind(jo_vec_sum$r, jo_vec_sum$eigen, jo_vec_sum$eigen_pval))
}
