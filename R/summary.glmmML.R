#' Summary of a glmmML object
#' 
#' It simply calls \code{print.glmmML}
#' 
#' 
#' @param object A glmmML object
#' @param \dots Additional arguments
#' @return Nothing is returned.
#' @note Preliminary
#' @author Göran Broström
#' @seealso \code{\link{print.glmmML}}
#' @keywords print
#' @export summary.glmmML
summary.glmmML <- function(object, ...){
    print.glmmML(object, ...)
}
