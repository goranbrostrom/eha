#' A summary of coxreg objects. 
#' 
#' 
#' 
#' 
#' @param object A \code{coxreg} object
#' @param \dots Additional ...
#' @author Göran Broström
#' @seealso \code{\link{print.coxreg}}
#' @keywords survival print
#' @examples
#' 
#' fit <- coxreg(Surv(enter, exit, event) ~ sex + civ, data = oldmort)
#' summary(fit)
#'  
#' @export
summary.coxreg <- function(object, ...){
    if (!object$nullModel){
        dr <- drop1(object, test = "Chisq")
        object$dr <- dr
        class(object) <- "summary.coxreg"
        coefficients <- cbind(object$coefficients, 
                              exp(object$coefficients),
                              sqrt(diag(object$var)))
        zval <- coefficients[, 1] / coefficients[, 3]
        pval <- pchisq(zval^2, df = 1, lower.tail = FALSE )
        coefficients <- cbind(coefficients, zval, pval)
        colnames(coefficients) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|")
        rownames(coefficients) <- names(object$coefficients)
    
    ##list(fit = object, coefficients = coefficients)
        object$coefficients <- coefficients
    }
    object
}
