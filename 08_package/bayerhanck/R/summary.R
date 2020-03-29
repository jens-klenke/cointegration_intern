#' Overview of the outcomes of
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
summary.co.test <- function(obj) {
  cat(c("----------------------------------------------------------",
        paste(c(obj$test, "Test"),
              collapse = " "),
        "----------------------------------------------------------"),
      sep = "\n")
  cat("Value of test statistic:" , obj$test.stat)
}

summary.bayerhanck <- function(obj) {
  cat(c("----------------------------------------------------------",
        "Bayerhanck Test for Non-Cointegration",
        "----------------------------------------------------------"),
      sep = "\n")
}
