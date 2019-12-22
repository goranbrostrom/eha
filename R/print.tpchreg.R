#' Prints tpchreg objects
#' 
#' More "pretty-printing"
#' 
#' Doesn't work with three-way or higher interactions.
#' 
#' @param x A \code{tpchreg} object, typically the result of running
#' \code{tpchreg}
#' @param digits Output format.
#' @param \dots Other arguments.
#' @return No value is returned.
#' @author Göran Broström
#' @seealso \code{\link{tpchreg}}, \code{\link{print.coxreg}}
#' @keywords survival
#' @export
print.tpchreg <-
 function(x, digits=max(options()$digits - 4, 3), ...){
    print.coxreg(x, digits, ...)
}
