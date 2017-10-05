#' The Lognormal Distribution
#' 
#' Density, distribution function, quantile function, hazard function,
#' cumulative hazard function, and random generation for the Lognormal distribution
#' with parameters \code{shape} and \code{scale}.
#' 
#' The Lognormal distribution with \code{scale} parameter \eqn{a} and \code{shape}
#' parameter \eqn{\sigma}{b} has hazard function given by \deqn{h(x) =
#' (b/\sigma)(x/\sigma)^(b-1)\exp((x / \sigma)^b)}{% h(x) =
#' (b/a)(x/a)^(b-1)exp((x / a)^b)} for \eqn{x \ge 0}{x >= 0}.
#' 
#' @name lnorm
#' @aliases lnorm dlnorm plnorm qlnorm hlnorm Hlnorm rlnorm
#' @usage 
#' hlnorm(x, meanlog = 0, sdlog = 1, shape = 1 / sdlog, scale = exp(meanlog),
#' prop = 1, log = FALSE) 
#' Hlnorm(x, meanlog = 0, sdlog = 1, shape = 1 / sdlog, scale = exp(meanlog), 
#' prop = 1, log.p = FALSE) 
#' @param x vector of quantiles.
#' @param meanlog mean in the Normal distribution.
#' @param sdlog,shape sdlog is standard deviation in the Normal distrimution, 
#' shape = 1/sdlog.
#' @param scale is exp(meanlog).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param prop proportionality constant in the extended Lognormal distribution.
#' @return \code{dlnorm} gives the density, \code{plnorm} gives the distribution
#' function, \code{qlnorm} gives the quantile function, \code{hlnorm} gives the
#' hazard function, \code{Hlnorm} gives the cumulative hazard function, and
#' \code{rlnorm} generates random deviates.
#' 
#' Invalid arguments will result in return value \code{NaN}, with a warning.
#' @keywords distribution
#' @export

hlnorm <- function(x, meanlog = 0, sdlog = 1,
                   shape = 1 / sdlog, scale = exp(meanlog), prop = 1,
                   log = FALSE){
    ## 'prop' added in version 2.3-0. 
    ## shape = 1 / sdlog, scale = exp(meanlog)
    meanlog <- log(scale)
    sdlog <- 1 / shape
    ret <- prop * dlnorm(x, meanlog, sdlog) /
        plnorm(x, meanlog, sdlog, lower.tail = FALSE)
    if (log) ret <- ifelse(ret <= 0, -Inf, log(ret))
    return (ret)

}

#' @export
Hlnorm <- function(x, meanlog = 0, sdlog = 1,
                   shape = 1 / sdlog, scale = exp(meanlog), prop = 1,
                   log.p = FALSE){
    ## 'prop' added in version 2.3-0. 
    ## shape = 1 / sdlog, scale = exp(meanlog)
    meanlog <- log(scale)
    sdlog <- 1 / shape
    ret <- -prop * plnorm(x, meanlog, sdlog, lower.tail = FALSE, log.p = TRUE)
    if (log.p) ret <- log(ret)
    return (ret)
}
