#' Prints a 'glmmML' object.
#' 
#' A glmmML object is the output of \code{glmmML}.
#' 
#' Nothing in particular.
#' 
#' @param x The glmmML object
#' @param digits Number of printed digits.
#' @param na.print How to print NAs
#' @param \dots Additional parameters, which are ignored.
#' @return A short summary of the object is printed.
#' @note This is the only summary method available for the moment.
#' @author Göran Broström
#' @seealso \code{\link{glmmML}}
#' @keywords print
#' @export print.glmmML
print.glmmML <- function(x,
                         digits = max(3, getOption("digits") - 3),
                         na.print = "",
                         ...){ 
    
    cat("\nCall: ", deparse(x$call), "\n\n")
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    coef <- x$coefficients
    se <- x$coef.sd

    tmp <- cbind(coef,
                 se,
                 coef/se,
                 signif(1 - pchisq((coef/se)^2, 1), digits - 1)
                 )
    dimnames(tmp) <- list(names(coef),
                          c("coef", "se(coef)", "z", "Pr(>|z|)")
                          )
    cat("\n")
    prmatrix(tmp)
    
    cat("\nScale parameter in mixing distribution: ", x$sigma,
        x$prior, "\n")
    cat("Std. Error:                             ", x$sigma.sd, "\n")

    pv <- 0.5 * pchisq(x$cluster.null.deviance - x$deviance, df = 1,
                       lower.tail = FALSE)
        cat("\n        LR p-value for H_0: sigma = 0: ", pv, "\n")
    if(x$boot){
        cat("\n Bootstrap p-value for H_0: sigma = 0: ",
        x$bootP, "(", x$boot, ")\n")
    }
    cat("\nResidual deviance:",
        format(signif(x$deviance, digits)), "on",
        x$df.residual, "degrees of freedom", 
        "\tAIC:",
        format(signif(x$aic, digits)), "\n")
}
