#' Prints summary.aftreg objects
#' 
#' @param x A \code{summary.aftreg} object, typically the result of running
#' \code{summary.aftreg}, summary on a phreg object.
#' @param digits Output format.
#' @param \dots Other arguments.
#' @return No value is returned.
#' @author Göran Broström
#' @seealso \code{\link{aftreg}}, \code{\link{summary.aftreg}}
#' @keywords survival
#' @export
print.summary.aftreg <- function(x, 
                                  digits = max(getOption("digits") - 3, 3),
                                  ...){
    class(x) <- c("summary.coxreg", class(x))
    print(x, digits)
}