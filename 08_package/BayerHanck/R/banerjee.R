#' Banerjee Test
#'
#' Executes Banerjee Test.
#'
#' @param formula An object of class "formula" to describe the model.
#' @param data An optional data frame containing the variables in the model.
#' @param lags Number of lags to be included.
#' @param trend Type of deterministic component to be inlcuded, "none" for no deterministics,
#' "const" for a constant and "trend" for a constant plus trend.
#'
#' @return \code{banerjee} returns an object of class \code{"co.test"}.
#' @export
#'
#' @references Banerjee, A., Dolado, J. J. and Mestre, R. (1998), Error-correction Mechanism Tests for Cointegration in a Single-equation Framework, Journal of Times Series Analysis 19(3), 267-283.
#'
#' @examples banerjee(linvestment ~ lincome + lconsumption, data = Lutkepohl)
banerjee <- function(formula, data, lags = 1, trend = "const"){

  #-----------------------------------------------------------------------------------------
  # Check Syntax
  #-----------------------------------------------------------------------------------------
  mf <- match.call()
  m <- match(c("formula", "data"), names(mf), 0L)
  if (is.null(data))
    stop()
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  mf <- na.omit(mf)
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)[, -1]
  trend <- match.arg(trend,
                     choices = c("none", "const", "trend"))

  #-----------------------------------------------------------------------------------------
  # Lag Matrix
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
  res <- matrix(NA, nrow = nrow(Xlag) - lags - 1, ncol = ncol(Xlag))

  if (identical(trend, "none")) {
    for (i in 1:ncol(Xlag)) {
      loop_lm <- lm(Hmisc::Lag(Xlag[, i], shift = 1)[-1] ~ W - 1)
      res[, i] <- as.numeric(loop_lm$residuals)
    }
  } else if (identical(trend, "const")) {
    for (i in 1:ncol(Xlag)) {
      loop_lm <- lm(Hmisc::Lag(Xlag[, i], shift = 1)[-1] ~ W)
      res[, i] <- as.numeric(loop_lm$residuals)
    }
  } else if (identical(trend, "trend")) {
    tr <- seq_along(y)[-1]
    W <- cbind(W, tr)
    for (i in 1:ncol(Xlag)) {
      loop_lm <- lm(Hmisc::Lag(Xlag[, i], shift = 1)[-1] ~ W)
      res[, i] <- as.numeric(loop_lm$residuals)
    }
    W <- W[, -which(colnames(W) == "tr")]
  }

  if (identical(trend, "none")) {
    BB_lm <- lm(Y_dif ~ W - 1)
  } else if (identical(trend, "const")) {
    BB_lm <- lm(Y_dif ~ W)
  } else if (identical(trend, "trend")) {
    tr <- seq_along(y)[-1]
    W <- cbind(W, tr)
    BB_lm <- lm(Y_dif ~ W)
    W <- W[, -which(colnames(W) == "tr")]
  }

  BB_res <- BB_lm$residuals

  if (identical(trend, "none")) {
    lm_res <- lm(BB_res ~ res -1)
  } else if (identical(trend, "const")) {
    lm_res <- lm(BB_res ~ res)
  } else if (identical(trend, "trend")) {
    tr <- seq_along(y)[-c(1, 2)]
    res <- cbind(res, tr)
    lm_res <- lm(BB_res ~ res)
    res <- res[, -which(colnames(res) == "tr")]
  }

  betas <- stats::coef(lm_res)
  var_mat <- vcov(lm_res)
  test.stat <- as.numeric(betas[2]/sqrt(var_mat[2, 2]))

  out <- list(test.stat = test.stat,
              lags = lags,
              trend = trend,
              betas = betas,
              var.cov = var_mat,
              test = "Banerjee",
              formula = formula)
  class(out) <- c("co.test", "list")
  cat(c("----------------------------------------------------------",
        "Banerjee Test",
        "----------------------------------------------------------",
        paste(c("Value of test statistic:", round(test.stat, 4)), collapse = " ")),
        sep = "\n")
  invisible(out)
}



