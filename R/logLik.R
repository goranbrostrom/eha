##
#        logLik methods
##

#' @export
logLik.coxreg <- function(object, ...){
    if (object$nullModel){
        out <- object$loglik[1]
        attr(out, "df") <- 0
    }else{
        out <- object$loglik[2]
        attr(out, "df") <- object$df
    }
    attr(out, "nobs") <- object$n.events
    class(out) <- "logLik"
    out
}

#' @export
logLik.phreg <- function(object, ...){
    out <- object$loglik[2]
    dd <- diag(object$var)
    attr(out, "df") <- sum(!is.na(dd) & dd > 0) #Stolen from logLik.survreg
    class(out) <- "logLik"
    out
}