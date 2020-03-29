#' Execute Banerjee Test
#'
#' @param formula An object of class "formula" to describe the model
#' @param data An optional data frame containing the variables in the model
#' @param lags Number of lags to be included
#' @param trend Type of deterministic component to be inlcuded, "none" for no deterministics,
#' "const" for a constant and "trend" for a constant plus trend
#'
#' @return
#' @export
#'
#' @examples
banerjee <- function(formula, data, lags = 1, trend = "const"){

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

  #-----------------------------------------------------------------------------------------
  # Lag Function
  #-----------------------------------------------------------------------------------------
  Xlag <- cbind(y, x)
  Y_dif <- diff(y) # muss als numeric vorliegen
  W <- diff(x)
  Xlag_diff <- as.data.frame(diff(Xlag))

  lag_mat <- matrix(NA, nrow = nrow(Xlag_diff), ncol = lags)
  lag_vec <- c()

  if (lags >= 1) {
    for (j in 1:ncol(Xlag_diff)) {
      for (i in 1:lags) {
        lag_mat[, i] <- Hmisc::Lag(Xlag_diff[, j], shift = i)
      }
      lag_vec <- cbind(lag_vec, lag_mat)
    }
    W <- cbind(W, lag_vec)
  }

  #-----------------------------------------------------------------------------------------
  # Banerjee Test
  #-----------------------------------------------------------------------------------------
  res_mat <- matrix(NA, nrow = nrow(Xlag) - lags - 1, ncol = ncol(Xlag))

  if (identical(trend, "none")) {
    for (i in 1:ncol(Xlag)) {
      loop_lm <- lm(Hmisc::Lag(Xlag[, i], shift = 1)[-1] ~ W - 1)
      res_mat[, i] <- as.numeric(loop_lm$residuals)
    }
  } else if (identical(trend, "const")) {
    for (i in 1:ncol(Xlag)) {
      loop_lm <- lm(Hmisc::Lag(Xlag[, i], shift = 1)[-1] ~ W)
      res_mat[, i] <- as.numeric(loop_lm$residuals)
    }
  } #else if (identical(trend, "trend")) {}

  if (identical(trend, "none")) {
    BB_lm <- lm(Y_dif ~ W - 1)
  } else if (identical(trend, "const")) {
    BB_lm <- lm(Y_dif ~ W)
  } #else if (identical(trend, "trend")) {}

  BB_res <- BB_lm$residuals

  if (identical(trend, "none")) {
    lm_res <- lm(BB_res ~ res_mat -1)
  } else if (identical(trend, "const")) {
    lm_res <- lm(BB_res ~ res_mat)
  } #else if (identical(trend, "trend)) {}

  betas <- stats::coef(lm_res)
  var_mat <- vcov(lm_res)
  test.stat <- as.numeric(betas[2]/sqrt(var_mat[2, 2]))
  names(test.stat) <- "banerjee"

  list(test.stat = test.stat)

  cat("hi")
}




