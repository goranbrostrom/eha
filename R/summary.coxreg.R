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
    dr <- drop1(object, test = "Chisq")
    object$dr <- dr
    class(object) <- "summary.coxreg"
    coefficients <- cbind(object$coefficients, 
                          sqrt(diag(object$var)))
    
    ##list(fit = object, coefficients = coefficients)
    object
}
