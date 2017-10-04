#' The EV Distribution
#' 
#' Density, distribution function, quantile function, hazard function,
#' cumulative hazard function, and random generation for the EV distribution
#' with parameters \code{shape} and \code{scale}.
#' 
#' The EV distribution with \code{scale} parameter \eqn{a} and \code{shape}
#' parameter \eqn{\sigma}{b} has hazard function given by \deqn{h(x) =
#' (b/\sigma)(x/\sigma)^(b-1)\exp((x / \sigma)^b)}{% h(x) =
#' (b/a)(x/a)^(b-1)exp((x / a)^b)} for \eqn{x \ge 0}{x >= 0}.
#' 
#' @name EV
#' @aliases EV dEV pEV qEV hEV HEV rEV
#' @usage dEV(x, shape = 1, scale = 1, log = FALSE) 
#' pEV(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) 
#' qEV(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) 
#' hEV(x, shape = 1, scale = 1, log = FALSE) 
#' HEV(x, shape = 1, scale = 1, log.p = FALSE) 
#' rEV(n, shape = 1, scale = 1)
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @param shape,scale shape and scale parameters, both defaulting to 1.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le
#' x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
#' @return \code{dEV} gives the density, \code{pEV} gives the distribution
#' function, \code{qEV} gives the quantile function, \code{hEV} gives the
#' hazard function, \code{HEV} gives the cumulative hazard function, and
#' \code{rEV} generates random deviates.
#' 
#' Invalid arguments will result in return value \code{NaN}, with a warning.
#' @keywords distribution
#' @export
pEV <- function(q, shape = 1, scale =  1,
                      lower.tail = TRUE, log.p = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    q <- ifelse(q <= 0, q, (q / scale)^shape)
    if (log.p){
        if (lower.tail){
            ret <- ifelse(q <= 0, -Inf, log1p(-exp(-expm1(q))))
        }else{
            ret <- ifelse(q <= 0, 0, -expm1(q))
        }
    }else{
        if (lower.tail){
            ret <- ifelse(q <= 0, 0, 1.0 - exp(-expm1(q)))
        }else{
            ret <- ifelse(q <= 0, 1, exp(-expm1(q)))
        }
    }

    return ( ret )
}

#' @export
dEV <- function(x, shape = 1, scale = 1, log = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    y <- (x / scale)^shape

    ret <- ifelse(x < 0, 0, (shape / scale) * (x / scale)^(shape - 1) *
                  exp(y - expm1(y)))

    if (log) ret <- log(ret)

    return ( ret )
}

#' @export
hEV <- function(x, shape = 1, scale = 1, log = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    y <- (x / scale)^shape

    ret <- ifelse(x < 0, 0,
                  (shape / scale) * (x / scale)^(shape - 1) * exp(y))

    if (log) ret <- log(ret)

    return ( ret )
}

#' @export
qEV <- function(p, shape = 1, scale = 1,
                      lower.tail = TRUE, log.p = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p

    ok <- (p >= 0) & (p <= 1)

    ret <- ifelse(ok, (1 / scale) * (log1p(-log1p(-p)))^(1 / shape), NaN)

    if (!all(ok)) warning("qEV produced NaN's")

    return ( ret )
}

#' @export
HEV <- function(x, shape = 1, scale = 1, log.p = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    y <- ifelse(x <= 0, 0, (x / scale)^shape)

    if (log.p){
        ret <- ifelse(y <= 0, -Inf, log(expm1(y)))
    }else{
        ret <- ifelse(y <= 0, 0, expm1(y))
    }

    return ( ret )
}

rEV <- function(n, shape = 1, scale = 1){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    y <- runif(n)

    return ( (1 / scale) * (log1p(-log1p(-y)))^(1 / shape) )
}
