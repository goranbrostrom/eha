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

#' @export
extractAIC.glmmML <- function(fit, scale, k = 2, ...) {
    if (k != 2) warning("Only k = 2 is implemented")
    edf <- length(fit$coefficients) + 1
    c(edf, fit$aic)
}

#' @export
extractAIC.glmmboot <- function(fit, scale, k = 2, ...) {
    if (k != 2) warning("Only k = 2 is implemented")
    edf <- length(fit$coefficients) + length(fit$frail)
    c(edf, fit$aic)
}

#' @export
nobs.glmmML <- function(object, ...){
    object$n
}

#' @export
nobs.glmmboot <- function(object, ...){
    object$n
}
