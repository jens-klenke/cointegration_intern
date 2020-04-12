#' @import urca
# https://github.com/cran/urca/blob/master/R/All-classes.R (Accessed: 11-09-2019)
setClass("ur.df.pvurt", representation(approxpval = "matrix"), contains="ur.df")

#' @import urca
# https://github.com/cran/urca/blob/master/R/All-classes.R (Accessed: 11-09-2019)
setClass("sumurca.pvurt", representation(approxpval = "otherornull"), contains = "sumurca")

# class for htest
setClass(
  "htest.pvurt",
  representation(
    statistic = "numeric",
    parameter = "numeric",
    alternative = "character",
    p.value = "numeric",
    pvurt = "numeric",
    method = "character",
    data.name = "character"
  )
)

