#' The Gompertz-Makeham Distribution
#' 
#' Density, distribution function, quantile function, hazard function,
#' cumulative hazard function, and random generation for the Gompertz-Makeham 
#' distribution with parameters \code{shape} and \code{scale}.
#' 
#' @details
#' The Gompertz-Makeham distribution with \code{shape} parameter \eqn{a} and \code{scale}
#' parameter \eqn{\sigma}{b} has hazard function given by 
#' \deqn{h(x) = a[2] + a[1] \exp(x/\sigma)}{h(x) = a[2] + a[1] exp(x/b)}
#' for \eqn{x \ge 0}{x >= 0}.
#' 
#' @name makeham
#' @aliases makeham dmakeham pmakeham qmakeham hmakeham Hmakeham rmakeham
#' @usage dmakeham(x, shape = c(1, 1), scale = 1, log = FALSE)
#' pmakeham(q, shape = c(1, 1), scale = 1, lower.tail = TRUE, log.p = FALSE)
#' qmakeham(p, shape = c(1, 1), scale = 1, lower.tail = TRUE, log.p = FALSE)
#' hmakeham(x, shape = c(1, 1), scale = 1, log = FALSE)
#' Hmakeham(x, shape = c(1, 1), scale = 1, log.p = FALSE)
#' rmakeham(n, shape = c(1, 1), scale = 1)
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @param shape A vector, default value c(1, 1).
#' @param scale  defaulting to 1.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le
#' x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
#' @return \code{dmakeham} gives the density, \code{pmakeham} gives the distribution
#' function, \code{qmakeham} gives the quantile function, \code{hmakeham} gives the
#' hazard function, \code{Hmakeham} gives the cumulative hazard function, and
#' \code{rmakeham} generates random deviates.
#' 
#' Invalid arguments will result in return value \code{NaN}, with a warning.
#' @keywords distribution
### The makeham-Makeham distribution, the hazard is
### h(t; shape, scale) = shape[2] + shape[1] * exp(t / scale), t >= 0
#' @export
pmakeham <- function(q, shape = c(1, 1), scale = 1,
                     lower.tail = TRUE, log.p = FALSE){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    if (all(q <= 0)){ # Fix trivial case first
        ret <- 1 - as.numeric(lower.tail)
        if (log.p) ret <- log(ret)
        return ( rep(ret, length(q)) )
    }
    ## Serious business:

    tmp <- ifelse(q <= 0,
                  0,
                  -( shape[2] * q + shape[1] * scale * expm1(q / scale) )
                  )
    if (lower.tail){
        if (log.p){
            ret <- ifelse(q <= 0, -Inf, log1p(-exp(tmp)))
        }else{
            ret <- ifelse(q <= 0, 0, -expm1(tmp))
        }
    }else{
        if (log.p){
            ret <- ifelse(q <= 0, 0, tmp)
        }else{
            ret <- ifelse(q <= 0, 1, exp(tmp))
        }
    }

    return ( ret )
}

#' @export
hmakeham <- function(x, shape = c(1, 1), scale = 1, log = FALSE){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    ret <- ifelse(x < 0,
                  0,
                  shape[2] + shape[1] * exp(x / scale))
    if (log) ret <- log(ret)

    return ( ret )
}

#' @export
dmakeham <- function(x, shape = c(1, 1), scale = 1, log = FALSE){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    ret <- ifelse(x < 0, 0, hmakeham(x, shape, scale, log = FALSE) *
                  pmakeham(x, shape, scale,
                           lower.tail = FALSE, log.p = FALSE))
    if (log) ret <- log(ret)

    return ( ret )
}

#' @export
Hmakeham <- function(x, shape = c(1, 1), scale = 1, log.p = FALSE){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    ret <- ifelse(x <= 0, 0,
                  shape[2] * x + shape[1] * scale * expm1(x / scale))
    if (log.p) ret <- log(ret)

    return ( ret )
}

#' @export
qmakeham <- function(p, shape = c(1, 1), scale = 1,
                     lower.tail = TRUE, log.p = FALSE){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    stop("Sorry, qmakeham is not yet implemented")
}

#' @export
rmakeham <- function(n, shape = c(1, 1), scale = 1){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    ## In the absence of 'qmakeham', we utilize the fact that a
    ## Gomperts-Makeham distributed random variable has the same
    ## distribution as the minimum of two independent random variables,
    ## one with an exponential and the other with a Gompertz distribution.

    x1 <- rexp(n, rate = shape[2])
    x2 <- rgompertz(n, shape = shape[1], scale = scale)

    return ( pmin(x1, x2) )
}
