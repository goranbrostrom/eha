#' Summary of a glmmboot object
#' 
#' It simply calls \code{print.glmmboot}
#' 
#' A summary method will be written soon.
#' 
#' @param object A glmmboot object
#' @param \dots Additional arguments
#' @return Nothing is returned.
#' @note Preliminary
#' @author Göran Broström
#' @seealso \code{\link{print.glmmboot}}
#' @keywords print
#' @export summary.glmmboot
summary.glmmboot <- function(object, ...){
    print.glmmboot(object, ...)
}
