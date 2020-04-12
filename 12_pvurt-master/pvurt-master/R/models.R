#' Models
#'
#' Models for approximating p-value.
#'
#' These are the functions based on polynomial regression and GAM to approximate
#' the p-value for ADF unit root tests of three types: without drift and trend
#' (woTD), with drift (wD) and with both drift and trend (wTD).
#'
#' @name models
#' @rdname models
#'
#' @param t Test statistics
#' @param n Sample size
#'
#' @return p-value from the ADF test statistic distributions.
NULL

#' @rdname models
woTD_poly3 <- function(t, n){

  woTD_t_min <- -5.042832
  woTD_t_max <- 4.679975

  t <- ifelse(t < woTD_t_min, woTD_t_min, t)
  t <- ifelse(t > woTD_t_max, woTD_t_max, t)

  p <- woTD_poly3_coef[1,1] +
    woTD_poly3_coef[2,1]*(t^1) +
    woTD_poly3_coef[3,1]*(t^2) +
    woTD_poly3_coef[4,1]*(t^3) +
    woTD_poly3_coef[5,1]*(1/n) +
    woTD_poly3_coef[6,1]*(t^1)*(1/n) +
    woTD_poly3_coef[7,1]*(t^2)*(1/n) +
    woTD_poly3_coef[8,1]*(t^3)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}

#' @rdname models
woTD_poly4 <- function(t, n){

  woTD_t_min <- -5.042832
  woTD_t_max <- 4.679975

  t <- ifelse(t < woTD_t_min, woTD_t_min, t)
  t <- ifelse(t > woTD_t_max, woTD_t_max, t)

  p <- woTD_poly4_coef[1,1] +
    woTD_poly4_coef[2,1]*(t^1) +
    woTD_poly4_coef[3,1]*(t^2) +
    woTD_poly4_coef[4,1]*(t^3) +
    woTD_poly4_coef[5,1]*(t^4) +
    woTD_poly4_coef[6,1]*(1/n) +
    woTD_poly4_coef[7,1]*(t^1)*(1/n) +
    woTD_poly4_coef[8,1]*(t^2)*(1/n) +
    woTD_poly4_coef[9,1]*(t^3)*(1/n) +
    woTD_poly4_coef[10,1]*(t^4)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}

#' @rdname models
woTD_poly5 <- function(t, n){

  woTD_t_min <- -5.042832
  woTD_t_max <- 4.679975

  t <- ifelse(t < woTD_t_min, woTD_t_min, t)
  t <- ifelse(t > woTD_t_max, woTD_t_max, t)

  p <- woTD_poly5_coef[1,1] +
    woTD_poly5_coef[2,1]*(t^1) +
    woTD_poly5_coef[3,1]*(t^2) +
    woTD_poly5_coef[4,1]*(t^3) +
    woTD_poly5_coef[5,1]*(t^4) +
    woTD_poly5_coef[6,1]*(t^5) +
    woTD_poly5_coef[7,1]*(1/n) +
    woTD_poly5_coef[8,1]*(t^1)*(1/n) +
    woTD_poly5_coef[9,1]*(t^2)*(1/n) +
    woTD_poly5_coef[10,1]*(t^3)*(1/n) +
    woTD_poly5_coef[11,1]*(t^4)*(1/n) +
    woTD_poly5_coef[12,1]*(t^5)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}

#' @rdname models
woTD_poly6 <- function(t, n){

  woTD_t_min <- -5.042832
  woTD_t_max <- 4.679975

  t <- ifelse(t < woTD_t_min, woTD_t_min, t)
  t <- ifelse(t > woTD_t_max, woTD_t_max, t)

  p <- woTD_poly6_coef[1,1] +
    woTD_poly6_coef[2,1]*(t^1) +
    woTD_poly6_coef[3,1]*(t^2) +
    woTD_poly6_coef[4,1]*(t^3) +
    woTD_poly6_coef[5,1]*(t^4) +
    woTD_poly6_coef[6,1]*(t^5) +
    woTD_poly6_coef[7,1]*(t^6) +
    woTD_poly6_coef[8,1]*(1/n) +
    woTD_poly6_coef[9,1]*(t^1)*(1/n) +
    woTD_poly6_coef[10,1]*(t^2)*(1/n) +
    woTD_poly6_coef[11,1]*(t^3)*(1/n) +
    woTD_poly6_coef[12,1]*(t^4)*(1/n) +
    woTD_poly6_coef[13,1]*(t^5)*(1/n) +
    woTD_poly6_coef[14,1]*(t^6)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}

#' @rdname models
#' @importFrom stats predict
#' @import mgcv
woTD_gam <- function(t, n){
  data <- data.frame(t = t, n = n)
  p <- predict(woTD_gam_model, newdata = data)
  p <- 1/(1 + exp(-p))
  return(p[[1]])
}

#' @rdname models
wD_poly3 <- function(t, n){

  wD_t_min <- -5.800194
  wD_t_max <- 3.0484

  t <- ifelse(t < wD_t_min, wD_t_min, t)
  t <- ifelse(t > wD_t_max, wD_t_max, t)

  p <- wD_poly3_coef[1,1] +
    wD_poly3_coef[2,1]*(t^1) +
    wD_poly3_coef[3,1]*(t^2) +
    wD_poly3_coef[4,1]*(t^3) +
    wD_poly3_coef[5,1]*(1/n) +
    wD_poly3_coef[6,1]*(t^1)*(1/n) +
    wD_poly3_coef[7,1]*(t^2)*(1/n) +
    wD_poly3_coef[8,1]*(t^3)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}

#' @rdname models
wD_poly4 <- function(t, n){

  wD_t_min <- -5.800194
  wD_t_max <- 3.0484

  t <- ifelse(t < wD_t_min, wD_t_min, t)
  t <- ifelse(t > wD_t_max, wD_t_max, t)

  p <- wD_poly4_coef[1,1] +
    wD_poly4_coef[2,1]*(t^1) +
    wD_poly4_coef[3,1]*(t^2) +
    wD_poly4_coef[4,1]*(t^3) +
    wD_poly4_coef[5,1]*(t^4) +
    wD_poly4_coef[6,1]*(1/n) +
    wD_poly4_coef[7,1]*(t^1)*(1/n) +
    wD_poly4_coef[8,1]*(t^2)*(1/n) +
    wD_poly4_coef[9,1]*(t^3)*(1/n) +
    wD_poly4_coef[10,1]*(t^4)*(1/n)
  p <- 1/(1+exp(-p))
  return(p[[1]])
}

#' @rdname models
wD_poly5 <- function(t, n){

  wD_t_min <- -5.800194
  wD_t_max <- 3.0484

  t <- ifelse(t < wD_t_min, wD_t_min, t)
  t <- ifelse(t > wD_t_max, wD_t_max, t)

  p <- wD_poly5_coef[1,1] +
    wD_poly5_coef[2,1]*(t^1) +
    wD_poly5_coef[3,1]*(t^2) +
    wD_poly5_coef[4,1]*(t^3) +
    wD_poly5_coef[5,1]*(t^4) +
    wD_poly5_coef[6,1]*(t^5) +
    wD_poly5_coef[7,1]*(1/n) +
    wD_poly5_coef[8,1]*(t^1)*(1/n) +
    wD_poly5_coef[9,1]*(t^2)*(1/n) +
    wD_poly5_coef[10,1]*(t^3)*(1/n) +
    wD_poly5_coef[11,1]*(t^4)*(1/n) +
    wD_poly5_coef[12,1]*(t^5)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}

#' @rdname models
wD_poly6 <- function(t, n){

  wD_t_min <- -5.800194
  wD_t_max <- 3.0484

  t <- ifelse(t < wD_t_min, wD_t_min, t)
  t <- ifelse(t > wD_t_max, wD_t_max, t)

  p <- wD_poly6_coef[1,1] +
    wD_poly6_coef[2,1]*(t^1) +
    wD_poly6_coef[3,1]*(t^2) +
    wD_poly6_coef[4,1]*(t^3) +
    wD_poly6_coef[5,1]*(t^4) +
    wD_poly6_coef[6,1]*(t^5) +
    wD_poly6_coef[7,1]*(t^6) +
    wD_poly6_coef[8,1]*(1/n) +
    wD_poly6_coef[9,1]*(t^1)*(1/n) +
    wD_poly6_coef[10,1]*(t^2)*(1/n) +
    wD_poly6_coef[11,1]*(t^3)*(1/n) +
    wD_poly6_coef[12,1]*(t^4)*(1/n) +
    wD_poly6_coef[13,1]*(t^5)*(1/n) +
    wD_poly6_coef[14,1]*(t^6)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}

#' @rdname models
#' @importFrom stats predict
#' @import mgcv
wD_gam <- function(t, n){
  data <- data.frame(t = t, n = n)
  p <- predict(wD_gam_model, newdata = data)
  p <- 1/(1 + exp(-p))
  return(p)
}


#' @rdname models
wTD_poly3 <- function(t, n){

  wTD_t_min <- -6.297081
  wTD_t_max <- 1.878922

  t <- ifelse(t < wTD_t_min, wTD_t_min, t)
  t <- ifelse(t > wTD_t_max, wTD_t_max, t)

  p <- wTD_poly3_coef[1,1] +
    wTD_poly3_coef[2,1]*(t^1) +
    wTD_poly3_coef[3,1]*(t^2) +
    wTD_poly3_coef[4,1]*(t^3) +
    wTD_poly3_coef[5,1]*(1/n) +
    wTD_poly3_coef[6,1]*(t^1)*(1/n) +
    wTD_poly3_coef[7,1]*(t^2)*(1/n) +
    wTD_poly3_coef[8,1]*(t^3)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}

#' @rdname models
wTD_poly4 <- function(t, n){

  wTD_t_min <- -6.297081
  wTD_t_max <- 1.878922

  t <- ifelse(t < wTD_t_min, wTD_t_min, t)
  t <- ifelse(t > wTD_t_max, wTD_t_max, t)

  p <- wTD_poly4_coef[1,1] +
    wTD_poly4_coef[2,1]*(t^1) +
    wTD_poly4_coef[3,1]*(t^2) +
    wTD_poly4_coef[4,1]*(t^3) +
    wTD_poly4_coef[5,1]*(t^4) +
    wTD_poly4_coef[6,1]*(1/n) +
    wTD_poly4_coef[7,1]*(t^1)*(1/n) +
    wTD_poly4_coef[8,1]*(t^2)*(1/n) +
    wTD_poly4_coef[9,1]*(t^3)*(1/n) +
    wTD_poly4_coef[10,1]*(t^4)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}

#' @rdname models
wTD_poly5 <- function(t, n){

  wTD_t_min <- -6.297081
  wTD_t_max <- 1.878922

  t <- ifelse(t < wTD_t_min, wTD_t_min, t)
  t <- ifelse(t > wTD_t_max, wTD_t_max, t)

  p <- wTD_poly5_coef[1,1] +
    wTD_poly5_coef[2,1]*(t^1) +
    wTD_poly5_coef[3,1]*(t^2) +
    wTD_poly5_coef[4,1]*(t^3) +
    wTD_poly5_coef[5,1]*(t^4) +
    wTD_poly5_coef[6,1]*(t^5) +
    wTD_poly5_coef[7,1]*(1/n) +
    wTD_poly5_coef[8,1]*(t^1)*(1/n) +
    wTD_poly5_coef[9,1]*(t^2)*(1/n) +
    wTD_poly5_coef[10,1]*(t^3)*(1/n) +
    wTD_poly5_coef[11,1]*(t^4)*(1/n) +
    wTD_poly5_coef[12,1]*(t^5)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}

#' @rdname models
wTD_poly6 <- function(t, n){

  wTD_t_min <- -6.297081
  wTD_t_max <- 1.878922

  t <- ifelse(t < wTD_t_min, wTD_t_min, t)
  t <- ifelse(t > wTD_t_max, wTD_t_max, t)

  p <- wTD_poly6_coef[1,1] +
    wTD_poly6_coef[2,1]*(t^1) +
    wTD_poly6_coef[3,1]*(t^2) +
    wTD_poly6_coef[4,1]*(t^3) +
    wTD_poly6_coef[5,1]*(t^4) +
    wTD_poly6_coef[6,1]*(t^5) +
    wTD_poly6_coef[7,1]*(t^6) +
    wTD_poly6_coef[8,1]*(1/n) +
    wTD_poly6_coef[9,1]*(t^1)*(1/n) +
    wTD_poly6_coef[10,1]*(t^2)*(1/n) +
    wTD_poly6_coef[11,1]*(t^3)*(1/n) +
    wTD_poly6_coef[12,1]*(t^4)*(1/n) +
    wTD_poly6_coef[13,1]*(t^5)*(1/n) +
    wTD_poly6_coef[14,1]*(t^6)*(1/n)
  p <- 1/(1+exp(-p))
  return(p)
}


#' @rdname models
#' @import mgcv
#' @importFrom stats predict
wTD_gam <- function(t, n){
  data <- data.frame(t = t, n = n)
  p <- predict(wTD_gam_model, newdata = data)
  p <- 1/(1 + exp(-p))
  return(p[[1]])
}
