#' The (Cumulative) Hazard Function of a Weibull Distribution
#' 
#' \code{hweibull} calculates the hazard function of a Weibull distribution,
#' and \code{Hweibull} calculates the corresponding cumulative hazard function.
#' 
#' See \link{dweibull}.
#' @name Weibull
#' @aliases hweibull Hweibull
#' @param x Vector of quantiles.
#' @param shape The shape parameter.
#' @param scale The scale parameter, defaults to 1.
#' @param log logical; if TRUE, the log of the hazard function is given.
#' @return The (cumulative) hazard function, evaluated at x.
#' @author Göran Broström
#' @seealso \code{\link{pweibull}}
#' @keywords survival
#' @examples
#' 
#' hweibull(3, 2, 1)
#' dweibull(3, 2, 1) / pweibull(3, 2, 1, lower.tail = FALSE)
#' Hweibull(3, 2, 1)
#' -pweibull(3, 2, 1, log.p = TRUE, lower.tail = FALSE)
#' 
#' @export
hweibull <- function(x, shape, scale = 1, log = FALSE){
    if (any(shape <= 0) || any(scale <= 0))
      stop("scale and shape must be positive")
    res <- ifelse(x < 0, 0, shape * (x / scale)^(shape - 1) / scale)
    if (log) res <- log(res)
    res
}

#' @export
Hweibull <- function(x, shape, scale = 1, log.p = FALSE){
    if (any(shape <= 0) || any(scale <= 0))
      stop("scale and shape must be positive")
    res <- ifelse(x < 0, 0, (x / scale)^shape)
    if (log.p) res <- log(res)
    res
}

