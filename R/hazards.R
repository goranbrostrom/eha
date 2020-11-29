#' Get baseline hazards atoms from a \code{coxreg} fit
#' 
#' @param x A \code{coxreg} object.
#' @return A list where each component is a two-column matrix representing
#' hazard atoms from one stratum. The first column is event time, and the 
#' second column is the corresponding hazard atom.
#' @export
hazards <- function(x){
    ## x is output from 'coxreg'.
    haz <- getHaz(x$y, strats = x$stratum, score = exp(x$linear.predictors))
    names(haz) <- x$strata
    haz
}