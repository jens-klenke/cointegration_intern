#' Johansen Test
#'
#' Executes Johansen Test.
#'
#' @param formula An object of class "formula" to describe the model.
#' @param data An optional data frame containing the variables in the model.
#' @param type Test to be conducted, either "eigen" or "trace".
#' @param lags Number of lags to be included.
#' @param trend Type of deterministic component to be inlcuded, "none" for no deterministics,
#' "const" for a constant and "trend" for a constant plus trend.
#'
#' @return \code{johansen} returns an object of class "co.test".
#' @export
#'
#' @examples
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
  mf <- na.omit(mf)
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

  out <- list(test.stat = test.stat,
              lags = lags,
              trend = trend,
              type = type,
              trace = jo_vec_sum[, 1:4],
              eigen = jo_vec_sum[, c(1, 5:6)])
  class(out) <- c("co.test", "list")
  out
}
