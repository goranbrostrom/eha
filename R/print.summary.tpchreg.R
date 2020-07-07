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
    class(x) <- c(class(x), "summary.coxreg")
    print.summary.coxreg(x, digits)
    cat("\nRestricted mean survival: ", x$rmean, "\n")
}