#' Prints phreg objects
#' 
#' This is the same as \code{\link{print.phreg}}
#' 
#' 
#' @param object A \code{phreg} object
#' @param \dots Additional ...
#' @author Göran Broström
#' @seealso \code{\link{print.coxreg}}
#' @keywords survival print
#' @examples
#' 
#' ## The function is currently defined as
#' function (object, ...) 
#' print(object)
#' 
#' @export
summary.phreg <- function(object, ...) print(object)
