#' Prints aftreg objects
#' 
#' 
#' @param object A \code{aftreg} object
#' @param \dots Additional ...
#' @author Göran Broström
#' @seealso \code{\link{print.coxreg}}
#' @keywords survival print
#' @examples
#' 
#' ## The function is currently defined as
#' function (object, ...) 
#' print(object)
#' 
#' @export
summary.aftreg <- function(object, ...){
    dr <- drop1(object, test = "Chisq")
    object$dr <- dr
    ncoef <- object$df
    ## Split coefficients into regression parameters and hazard ditto.
    if (!(object$dist == "pch")){
        hazards <- object$coefficients[-(1:ncoef)]
    }
    coefficients <- object$coefficients[1:ncoef]
    
    ## Regression parameters:
    rawnames <- names(coefficients)
    varcoef <- diag(object$var[1:ncoef, 1:ncoef, drop = FALSE])
    varhaz <- diag(object$var[-(1:ncoef), -(1:ncoef), drop = FALSE])
    class(object) <- "summary.aftreg"
    coefficients <- cbind(coefficients, 
                          exp(coefficients),
                          sqrt(varcoef))
    zval <- coefficients[, 1] / coefficients[, 3]
    pval <- pchisq(zval^2, df = 1, lower.tail = FALSE )
    coefficients <- cbind(coefficients, zval, pval)
    colnames(coefficients) <- c("coef", "exp(coef)", "se(coef)", "z", "Wald p"
    )
    rownames(coefficients) <- rawnames
    
    ## Hazard parameters:
    if (!(object$dist == "pch")){
        haznames <- names(hazards)
        hazards <- cbind(hazards, sqrt(varhaz))
        colnames(hazards) <- c("par", "se(par)")
        rownames(hazards) <- haznames
    }
    ##list(fit = object, coefficients = coefficients)
    object$coefficients <- coefficients
    object$hazards <- hazards
    object
    
}
