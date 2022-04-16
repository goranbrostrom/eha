#' A summary of phreg objects. 
#' 
#' 
#' 
#' 
#' @param object A \code{phreg} object
#' @param \dots Additional ...
#' @author Göran Broström
#' @seealso \code{\link{print.phreg}}
#' @keywords survival 
#' @examples
#' 
#' fit <- phreg(Surv(enter, exit, event) ~ sex + civ, 
#' data = oldmort[oldmort$region == "town", ])
#' summary(fit)
#'  
#' @export
summary.phreg <- function(object, ...){
    if (!object$nullModel){
        ncoef <- object$df
        dr <- drop1(object, test = "Chisq")
        object$dr <- dr
        ## Split coefficients into regression parameters and hazard ditto.
        if (!(object$dist == "pch")){
            hazards <- object$coefficients[-(1:ncoef)]
        }
        coefficients <- object$coefficients[1:ncoef]
    }else{ # null Model
        hazards <- object$coefficients
    }

    
    ## Regression parameters:
    if (!object$nullModel){
        rawnames <- names(coefficients)
        varcoef <- diag(object$var[1:ncoef, 1:ncoef, drop = FALSE])
        varhaz <- diag(object$var[-(1:ncoef), -(1:ncoef), drop = FALSE])

        coefficients <- cbind(coefficients, 
                              exp(coefficients),
                              sqrt(varcoef))
        zval <- coefficients[, 1] / coefficients[, 3]
        pval <- pchisq(zval^2, df = 1, lower.tail = FALSE )
        coefficients <- cbind(coefficients, zval, pval)
        colnames(coefficients) <- c("coef", "exp(coef)", "se(coef)", "z", "Wald p")
        rownames(coefficients) <- rawnames
        object$coefficients <- coefficients
    }
    ## Hazard parameters:
    if (!(object$dist == "pch")){
        haznames <- names(hazards)
        hazards <- cbind(hazards, sqrt(varhaz))
        colnames(hazards) <- c("par", "se(par)")
        rownames(hazards) <- haznames

        object$hazards <- hazards
    }
    
    class(object) <- "summary.phreg"

    ##list(fit = object, coefficients = coefficients)
    object
}
