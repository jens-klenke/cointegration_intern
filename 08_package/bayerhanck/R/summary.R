summary.co.test <- function(object) {
  cat(c("----------------------------------------------------------",
        paste(c(object$test, "Test"),
              collapse = " "),
        "----------------------------------------------------------",
        paste(c("Formula:", object$formula),
              collapse = " "),
        paste(c("Lags:", object$lags),
              collapse = " "),
        paste(c("Trend:", object$trend),
              collapse = " "),
        " ",
        paste(c("Value of test statistic:", round(object$test.stat, 4)),
              collapse = " ")),
      sep = "\n")
  out <- list(bh.test = object$bh.test)
  class(out) <- c("summary.co.test")
  invisible(out)
}


summary.bh.test <- function(object) {
  cat(c("----------------------------------------------------------",
      "Bayerhanck Test for Non-Cointegration",
      "----------------------------------------------------------",
      paste(c("Formula:", object$formula),
            collapse = " "),
      paste(c("Lags:", object$lags),
            collapse = " "),
      paste(c("Trend:", object$trend),
            collapse = " "),
      " ",
      "Underlying Tests:"),
    sep = "\n")
  test.mat <- as_tibble(rbind(object$test.stat,
                                  object$pval.stat))
  rownames(test.mat) <- c("Test Statistics", "p-Values")
  print(test.mat)
  cat(c(" ",
        paste(c("Value of the Fisher Type Test statistic:", round(object$bh.test, 4)),
              collapse = " ")),
      sep = "\n")
  crit <- if (identical(object$crit, 0.01) | identical(object$crit, 0.05)) {
    substr(object$crit, 4, 4)
  } else {
    paste(c(substr(object$crit, 3, 3), 0),
          collapse = "")
  }
  cat(c(paste(
    paste(c(crit, "%"), collapse = ""),
    "critical value for Test statistic:", object$crit.val)))
  out <- list(bh.test = object$bh.test)
  class(out) <- c("summary.bh.test")
  invisible(out)
}




