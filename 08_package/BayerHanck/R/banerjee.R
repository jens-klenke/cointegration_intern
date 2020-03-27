banerjee <- function(formula, data, lags){

  mf <- match.call()
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  offset <- model.offset(mf)
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)[, -1]

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

  ### Loop
  res_mat <- matrix(NA, nrow = nrow(Xlag) - lags - 1, ncol = ncol(Xlag))

  for (i in 1:ncol(Xlag)) {
    loop_lm <- lm(Hmisc::Lag(Xlag[, i], shift = 1)[-1] ~ W)
    res_mat[, i] <- as.numeric(loop_lm$residuals)
  }

  ### Boswijk Test
  BB_lm <- lm(Y_dif ~ W)
  BB_res <- BB_lm$residuals

  lm_res <- lm(BB_res ~ res_mat)
  betas <- coef(lm_res)
  var_mat <- vcov(lm_res)

  # Banerjee Test
  test.stat <- as.numeric(betas[2]/sqrt(var_mat[2, 2]))
  print(test.stat)
}

banerjee(data = df, formula = linvestment ~ lincome + lconsumption, lags = 1)
