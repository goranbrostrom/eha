#' The Gompertz Distribution
#' 
#' Density, distribution function, quantile function, hazard function,
#' cumulative hazard function, and random generation for the Gompertz 
#' distribution with parameters \code{shape} and \code{scale}.
#' 
#' @details
#' The Gompertz distribution with \code{scale} parameter \eqn{a} and \code{shape}
#' parameter \eqn{\sigma}{b} has hazard function given by 
#' \deqn{h(x) = a \exp(x/\sigma)}{%
#' h(x) = a exp(x/b)}
#' for \eqn{x \ge 0}{x >= 0}.
#' If \code{param = "canonical"}, then then a --> a/b, so that b is a
#' true scale parameter (for any fixed a), and b is an 'AFT parameter'.
#' If \code{param = "rate"}, then b --> 1/b.
#' @name Gompertz
#' @aliases gompertz dgompertz pgompertz qgompertz hgompertz Hgompertz rgompertz
#' @usage dgompertz(x, shape = 1, scale = 1, rate, log = FALSE, 
#' param = c("default", "canonical", "rate")) 
#' pgompertz(q, shape = 1, scale = 1, rate, lower.tail = TRUE, log.p = FALSE, 
#' param = c("default", "canonical", "rate")) 
#' qgompertz(p, shape = 1, scale = 1, rate, lower.tail = TRUE, log.p = FALSE, 
#' param = c("default", "canonical", "rate")) 
#' hgompertz(x, shape = 1, scale = 1, rate, log = FALSE, 
#' param = c("default", "canonical", "rate")) 
#' Hgompertz(x, shape = 1, scale = 1, rate, log.p = FALSE, 
#' param = c("default", "canonical", "rate")) 
#' rgompertz(n, shape = 1, scale = 1, , rate, 
#' param = c("default", "canonical", "rate"))
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @param shape,scale shape and scale parameters, both defaulting to 1.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le
#' x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
#' @param param default or canonical or rate.
#' @return \code{dgompertz} gives the density, \code{pgompertz} gives the distribution
#' function, \code{qgompertz} gives the quantile function, \code{hgompertz} gives the
#' hazard function, \code{Hgompertz} gives the cumulative hazard function, and
#' \code{rgompertz} generates random deviates.
#' 
#' Invalid arguments will result in return value \code{NaN}, with a warning.
#' @keywords distribution
#' @export
pgompertz <- function(q, shape = 1, scale = 1, rate,
                     lower.tail = TRUE, log.p = FALSE,
                      param = c("default", "canonical", "rate")){

    n <- length(q)
    if (any(scale == 0)){
        if (log.p) return(rep(-Inf, n))
        else return(rep(0, n))
    }
    ##if ( any(c(shape, scale) <= 0) ){
    if ( any(scale < 0) ){
        cat("scale = ", scale, "\n")
        warning("Negative scale")
        return(NaN)
    }

    if ( any(shape < 0) ){
        warning("Negative shape")
        return(NaN)
    }

    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param == "rate") {
        scale <- 1 / rate
    }else if (param != "default") stop("Illegal 'param'")
    ##

    y <- ifelse(q <= 0, 0, -shape * scale * expm1(q / scale))

    if (log.p){
        if (lower.tail){
            ret <- ifelse(q <= 0,
                          -Inf,
                          log(-expm1(y))
                          )
        }else{
            ret <- ifelse(q <= 0,
                          0,
                          y
                          )
        }
    }else{
        if (lower.tail){
            ret <- ifelse(q <= 0,
                          0,
                          -expm1(y) # = 1 - exp(y)
                          )
        }else{
            ret <- ifelse(q <= 0,
                          1,
                          exp(y)
                          )
        }
    }

    return ( ret )
}

#' @export
dgompertz <- function(x, shape = 1, scale = 1, rate, log = FALSE,
                      param = c("default", "canonical", "rate")){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    
    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param == "rate"){
        scale = 1 / rate
    }else if (param != "default") stop("Illegal 'param'")
    ##

    y <- ifelse(x < 0, 0, x / scale)

    if (log){
        ret <- ifelse(x < 0,
                      0,
                      log(shape) + y - shape * scale * expm1(y))
    }else{
        ret <- ifelse(x < 0,
                      0,
    ##                  shape * exp(y) * exp(-shape * scale * expm1(y)))
                      shape * exp(y - shape * scale * expm1(y)))
    }

    return ( ret )
}

#' @export
hgompertz <- function(x, shape = 1, scale = 1, rate, log = FALSE,
                      param = c("default", "canonical", "rate")){

    if ( any(shape < 0) ){
        cat("shape = ", shape, "\n")
        stop("Negative shape")
    }
    
    param <- param[1]

    if (param == "rate"){
        if (missing(rate)) rate <- 1 / scale
        ret <- shape * exp(rate * x)
        if (log) ret <- log(ret)
        return(ret)
    }
    
    ##if ( any(c(shape, scale) <= 0) ){
    if ( any(scale <= 0) ){
        cat("scale = ", scale, "\n")
        warning("Non-positive scale")
        return(rep(Inf, length(x)))
    }


    ## New in 2.1-3:

    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param != "default") stop("Illegal 'param'")
    ##

    
    if (log) {
        ret <- ifelse(x < 0,
                      -Inf,
                      log(shape) + x / scale
                      )
    }else{
        ret <- ifelse(x < 0,
                      0,
                      shape * exp(x / scale)
                      )
    }

    return ( ret )
}

#' @export
qgompertz <- function(p, shape = 1, scale = 1, rate,
                     lower.tail = TRUE, log.p = FALSE,
                      param = c("default", "canonical", "rate")){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param == "rate"){
        scale <- 1 / rate
    }else if (param != "default") stop("Illegal 'param'")
    ##
    
    if (log.p) p <- exp(p)

    ok <- (p >= 0) & (p <= 1)


    ret <- ifelse(ok,
                  scale * log1p(-log1p(-p) / (shape * scale)),
                  NaN)

    if (!all(ok)) warning("qgompertz produced NaN's")

    return ( ret )
}

#' @export
Hgompertz <- function(x, shape = 1, scale = 1, rate, log.p = FALSE,
                      param = c("default", "canonical", "rate")){

    if ( any(shape < 0) ){
        cat("shape = ", shape, "\n")
        stop("Negative shape")
    }
    
    param <- param[1]
    
    if (param == "rate"){
        if (missing(rate)) rate <- 1 / scale
        if (isTRUE(all.equal(rate, 0))){
            ret <- shape * x
        }else{
            y <- x * rate
            ret <- ifelse(x <= 0,
                          0,
                          shape / rate * expm1(y)
            )
        }
        if (log.p) ret <- log(ret)
        return(ret)
    }
    
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param != "default") stop("Illegal 'param'")
    ##

    
    y <- x / scale

    ret <- ifelse(x <= 0,
                  0,
                  shape * scale * expm1(y)
                  )
    if (log.p) ret <- log(ret)

    return ( ret )
}

#' @export
rgompertz <- function(n, shape = 1, scale = 1, rate,
                      param = c("default", "canonical", "rate")){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param == "rate"){
        scale <- 1 / rate
    }else if (param != "default") stop("Illegal 'param'")
    ##

    
    y <- runif(n)

    return ( qgompertz(y, shape, scale, param = "default") )# Note!
}
