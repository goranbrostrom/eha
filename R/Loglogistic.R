#' The Loglogistic Distribution
#' 
#' Density, distribution function, quantile function, hazard function,
#' cumulative hazard function, and random generation for the Loglogistic distribution
#' with parameters \code{shape} and \code{scale}.
#' 
#' The Loglogistic distribution with \code{scale} parameter \eqn{a} and \code{shape}
#' parameter \eqn{\sigma}{b} has hazard function given by \deqn{h(x) =
#' (b/\sigma)(x/\sigma)^(b-1)\exp((x / \sigma)^b)}{% h(x) =
#' (b/a)(x/a)^(b-1)exp((x / a)^b)} for \eqn{x \ge 0}{x >= 0}.
#' 
#' @name Loglogistic
#' @aliases Llogis dllogis pllogis qllogis hllogis Hllogis rllogis
#' @usage dllogis(x, shape = 1, scale = 1, log = FALSE) 
#' pllogis(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) 
#' qllogis(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) 
#' hllogis(x, shape = 1, scale = 1, prop = 1, log = FALSE) 
#' Hllogis(x, shape = 1, scale = 1, prop = 1, log.p = FALSE) 
#' rllogis(n, shape = 1, scale = 1)
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @param shape,scale shape and scale parameters, both defaulting to 1.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le
#' x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
#' @param prop proportionality constant in the extended Loglogistic distribution.
#' @return \code{dllogis} gives the density, \code{pllogis} gives the distribution
#' function, \code{qllogis} gives the quantile function, \code{hllogis} gives the
#' hazard function, \code{Hllogis} gives the cumulative hazard function, and
#' \code{rllogis} generates random deviates.
#' 
#' Invalid arguments will result in return value \code{NaN}, with a warning.
#' @keywords distribution
#' @export
pllogis <- function(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    location <- log(scale)
    scale <- 1 / shape
    if (log.p){
        if (lower.tail){
            ret <- ifelse(q <= 0,
                          -Inf,
                          plogis(log(q), location, scale, lower.tail, log.p)
                          )
        }else{
            ret <- ifelse(q <= 0,
                          0,
                          plogis(log(q), location, scale, lower.tail, log.p)
                          )
        }
    }else{ #! log.p
        if (lower.tail){
            ret <- ifelse(q <= 0,
                          0,
                          plogis(log(q), location, scale, lower.tail, log.p)
                          )
        }else{
            ret <- ifelse(q <= 0,
                          1,
                          plogis(log(q), location, scale, lower.tail, log.p)
                          )
        }
    }

    return ( ret )
}

#' @export
dllogis <- function(x, shape = 1, scale = 1, log = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }
    y <- x / scale


    ret <- ifelse(y < 0,
                  0,
                  (shape / scale) * y^(shape - 1) / (1 + y^shape)^2
                  )
    if (log) ret <- ifelse(ret <= 0, -Inf, log(ret))

    return (ret)
}

#' @export
hllogis <- function(x, shape = 1, scale = 1, prop = 1, log = FALSE){
    ## New argument 'prop', a proportionality parameter in
    ## an extended loglogistic distribution, added in version 2.3-0.
    if ( any(c(shape, scale, prop) <= 0) ){
        warning("Non-positive shape, scale, or prop")
        return(NaN)
    }
    y <- x / scale
    ret <- prop * ifelse(x < 0,
                         0,
                         (shape / scale) * y^(shape - 1) / (1 + y^shape)
                         )
    if (log) ret <- ifelse(ret <= 0, -Inf, log(ret))

    return(ret)
}

#' @export
qllogis <- function(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p

    location <- log(scale)
    scale <- 1 / shape

    ## Let 'qlogis do the checking....
    ret <- qlogis(p, location, scale, lower.tail, log.p)

    return(exp(ret))
}


#' @export
Hllogis <- function(x, shape = 1, scale = 1, prop = 1, log.p = FALSE){
    ## New argument 'prop', a proportionality parameter in
    ## an extended loglogistic distribution, added in version 2.3-0.
    if ( any(c(shape, scale, prop) <= 0) ){
        warning("Non-positive shape, scale, or prop")
        return(NaN)
    }

    ret = -prop * pllogis(x, shape, scale, lower.tail = FALSE, log.p = TRUE)
    if (log.p) ret <- log(ret)

    return (ret)
}

#' @export
rllogis <- function(n, shape = 1, scale = 1){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }
    location <- log(scale)
    scale <- 1 / shape

    exp(rlogis(n, location, scale))
}
