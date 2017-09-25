#' Prints coxreg objects
#' 
#' This is the same as \code{\link{print.coxreg}}
#' 
#' 
#' @param object A \code{coxreg} object
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
#' @export summary.coxreg
summary.coxreg <- function(object, ...) print(object)
