#' Prints summary.phreg objects
#' 
#' @param x A \code{summary.phreg} object, typically the result of running
#' \code{summary.phreg}, summary on a phreg object.
#' @param digits Output format.
#' @param \dots Other arguments.
#' @return No value is returned.
#' @author Göran Broström
#' @seealso \code{\link{phreg}}, \code{\link{summary.phreg}}
#' @keywords survival
#' @export
print.summary.phreg <- function(x, 
                                  digits = max(getOption("digits") - 3, 3),
                                  ...){
    class(x) <- c("summary.coxreg", class(x))
    print(x, digits)
}