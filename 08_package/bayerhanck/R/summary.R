summary.co.test <- function(object) {
  cat(c("----------------------------------------------------------",
        paste(c(object$test, "Test"),
              collapse = " "),
        "----------------------------------------------------------"),
      sep = "\n")
  cat("Value of test statistic:" , obj$test.stat)
}

summary.bayerhanck <- function(object) {
cat(c("----------------------------------------------------------",
      "Bayerhanck Test for Non-Cointegration",
      "----------------------------------------------------------"),
    sep = "\n")
}
