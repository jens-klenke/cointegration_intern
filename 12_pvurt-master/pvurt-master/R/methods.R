# https://github.com/cran/urca/blob/master/R/methods-show.R (Accessed: 11-09-2019)
show.urca.pvurt <- function(object){
  title <- paste("#", object@test.name, "Unit Root / Cointegration Test #", sep=" ")
  row <- paste(rep("#", nchar(title)), collapse="")
  ptext <- paste0("Approximated P value is (pvurt):", paste(rep(" ", 3), collapse = ""))
  cat("\n")
  cat(row, "\n")
  cat(title, "\n")
  cat(row, "\n")
  cat("\n")
  cat("The value of the test statistic is:", round(object@teststat, 4), "\n")
  cat(ptext, round(object@approxpval, 4), "\n")
  cat('\n')
}

# https://github.com/cran/urca/blob/master/R/methods-show.R (Accessed: 11-09-2019)
setMethod("show", "ur.df.pvurt", show.urca.pvurt)

# https://github.com/cran/urca/blob/master/R/methods-summary.R (Accessed: 11-09-2019)
setMethod("summary", "ur.df.pvurt", function(object){
  return(new("sumurca.pvurt", classname="ur.df.pvurt", test.name=object@test.name,
             testreg=object@testreg, teststat=object@teststat, approxpval=object@approxpval,
             cval=object@cval, bpoint=NULL, signif=NULL, model=object@model,
             type=NULL, auxstat=NULL, lag=NULL, H=NULL, A=NULL, lambda=NULL,
             pval=NULL, V=NULL, W=NULL, P=NULL))
})

# https://github.com/cran/urca/blob/master/R/methods-show.R (Accessed: 11-09-2019)
setMethod("show", "sumurca.pvurt", function(object){
  if(object@classname=="ur.df.pvurt"){
    title <- paste("#", object@test.name, "Unit Root Test #", sep=" ")
    row <- paste(rep("#", nchar(title)), collapse="")
    ptext <- "Approximated P value is (pvurt):"
    cat("\n")
    cat(row, "\n")
    cat(title, "\n")
    cat(row, "\n")
    cat("\n")
    cat('Test regression', object@model, '\n')
    cat('\n')
    print(object@testreg)
    cat('\n')
    cat('Value of test-statistic is:    ', round(object@teststat, 4), '\n')
    cat(ptext, round(object@approxpval, 4), "\n")
    cat('\n')
    cat('Critical values for test statistics: \n')
    print(object@cval)
    cat('\n')
    invisible(object)
  }
})


show.htest.pvurt <- function(object){
  digits <- 4
  cat("\n")
  cat(strwrap(object@method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data:  ", object@data.name, "\n", sep = "")
  out <- character()
  if (!is.null(object@statistic))
    out <- c(out, paste(names(object@statistic), "=", format(signif(object@statistic,
                                                                    max(1L, digits)))))
  if (!is.null(object@parameter))
    out <- c(out, paste(names(object@parameter), "=", format(signif(object@parameter,
                                                                    max(1L, digits)))))
  if (!is.null(object@p.value)) {
    fp <- format.pval(object@p.value, digits = max(1L, digits))
    out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("=", fp)))
  }
  if (!is.null(object@pvurt)) {
    fpp <- format.pval(object@pvurt, digits = max(1L, digits))
    out <- c(out, paste("pvurt", if (substr(fpp, 1L, 1L) == "<") fpp else paste("=", fpp)))
  }

  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("alternative hypothesis: ")
  cat(object@alternative, "\n", sep = "")

  cat("\n")

}

setMethod("show", "htest.pvurt", show.htest.pvurt)

