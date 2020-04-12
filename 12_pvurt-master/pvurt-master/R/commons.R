#' Common functions
#'
#' Common functions used by top level functions.
#'
#' \code{modelName} generates the name of function based on the arguments
#' \code{model}, \code{type} and \code{d} (when \code{model = "poly"}).\cr
#' \code{computePValueMain} is called by the other top level functions to
#' approximate the p-value.
#'
#' @name commons
#'
#' @rdname commons
#'
#' @param object Numeric value or an object (\code{fHTEST}, \code{ur.df} or
#'   \code{htest}) for which p-value needs to
#'   approximated.
#' @param n Sample size.
#' @param type The type of unit root test. Currently supports: \code{nc} for
#'   test without drift and trend, \code{c} for test with only drift and
#'   \code{ct} for test with both drift and trend.
#' @param model The model type to be used for approximation. Available is GAM
#'   and polynomial regression. If \code{gam} is chosen, then \code{d} has no
#'   effect.
#' @param d The degree for polynomial. \code{d} must be \eqn{\geq 3} and
#'   \eqn{\leq 6}. If \code{gam} is chosen, then \code{d} has no effect.
NULL

#' @rdname commons
computePValueMain <- function(
  object, n, model = c("gam", "poly"), type = c("nc", "c", "ct"), d = NULL
){
  model <- match.arg(model)
  type <- match.arg(type)
  t <- object
  n <- n
  if(model == "poly"){
    if(is.null(d)){
      stop("Specify a value for d when 'poly' is selected.")
    } else if (!is.null(d) & (d < 3 | d > 6)){
      stop("d must be >= 3 and  =< 6.")
    }else {
      d <- d
    }
  }

  if(type == "nc" & model == "poly" & n < 50){
    warning("It is recommended to use the GAM model instead, as the sample size can be too small for polynomial model.")
  }

  model_ <- modelName(type = type, model = model, d = d)
  p_ <- do.call(model_, list(t=t, n=n))

  return(p_)
}

#' @rdname commons
modelName <- function(type, model, d){
  typeVector <- c("woTD", "wD", "wTD")
  names(typeVector) <- c("nc", "c", "ct")

  modelName <- paste0(typeVector[[type]], "_", model)
  if(model == "poly") modelName <- paste0(modelName, d)

  return(modelName)
}
