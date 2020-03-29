#' Summary method for objects from bayerhanck
#'
#' @usage
#' ## S3 method for class "co.test"
#' ## S3 method for class "bayerhanck"
#'
#' @param object An object of class \code{co.test} or \code{bayerhanck}.
#'
#' @return
#' @export
#'
#' @examples
summary.co.test <- function(object) {
  cat(c("----------------------------------------------------------",
        paste(c(object$test, "Test"),
              collapse = " "),
        "----------------------------------------------------------"),
      sep = "\n")
  cat("Value of test statistic:" , obj$test.stat)
}
