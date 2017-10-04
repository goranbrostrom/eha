#' @importFrom stats nobs extractAIC
#' @export 
extractAIC.coxreg <- function(fit, scale, k = 2, ...) 
{
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

#' @export
extractAIC.phreg <- function(fit, scale, k = 2, ...) 
{
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

#' @export
extractAIC.aftreg <- function(fit, scale, k = 2, ...) 
{
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

#' @export
nobs.coxreg <- function(object, ...){
    object$n
}

#' @export
nobs.phreg <- function(object, ...){
    object$n
}

#' @export
nobs.aftreg <- function(object, ...){
    object$n
}
