#' Prints aftreg objects
#' 
#' This is the same as \code{\link{print.aftreg}}
#' 
#' 
#' @param object A \code{aftreg} object
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
#' @export summary.aftreg
summary.aftreg <- function(object, ...) print(object)
