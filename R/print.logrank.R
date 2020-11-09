#' Prints logrank objects
#' 
#' The result of \code{logrank} is printed
#' 
#' 
#' @param x A \code{logrank} object
#' @param digits Precision in printing
#' @param \dots Not used.
#' @return The input is returned invisibly.
#' @author Göran Broström
#' @seealso \code{\link{logrank}}, \code{\link{coxreg}}
#' @keywords survival regression
#' @export
print.logrank <- function(x, digits=max(options()$digits - 4, 6), ...){
  
    if (class(x) != "logrank") stop("Not a 'logrank' object.")
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    cat("\n     The log-rank test\n") 
    if (!is.null(cl <- x$call)){
        cat("\nCall:\n")
        dput(cl)
        cat("\n")
    }
    cat(" X-squared = ", x$test.statistic, ", df = ", x$df,
        ", p-value = ", format.pval(x$p.value, digits = digits), "\n")
    invisible(x)
}
