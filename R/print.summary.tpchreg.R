#' Prints summary.tpchreg objects
#' 
#' @param x A \code{summary.tpchreg} object, typically the result of running
#' \code{summary.tpchreg}, summary on a tpchreg object.
#' @param digits Output format.
#' @param \dots Other arguments.
#' @return No value is returned.
#' @author Göran Broström
#' @seealso \code{\link{tpchreg}}, \code{\link{summary.tpchreg}}
#' @keywords survival
#' @export
print.summary.tpchreg <- function(x, 
                                  digits = max(getOption("digits") - 3, 3),
                                  ...){
    class(x) <- c("summary.coxreg", class(x))
    print(x, digits)
    ivl <- paste("(", min(x$cuts), ", ", max(x$cuts), "]", sep = "")
    if (x$n.strata == 1){
        cat("\nRestricted mean survival: ", x$rmean, "in", ivl, "\n")
    }else{
        names(x$rmean) <- x$strata
        cat("\nRestricted mean survival in", ivl, ": \n")
        print(x$rmean)
        cat("\n")
    }
}