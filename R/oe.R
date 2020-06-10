#' Create an oe object
#'
#' Create an \emph{oe} ("occurrence/exposure") object, used as a response 
#' variable in a model formula specifically in \code{\link{tpchreg}}. 
#'
#' @usage oe(count, exposure)
#' @param count Number of events, a non-negative integer-valued vector.
#' @param exposure exposure time corresponding to count. 
#' A positive numeric vector.
#'
#'@seealso \code{\link{tpchreg}}.
#'@export
oe <- function(count, exposure){
    if (!is.numeric(count)) stop("count not numeric.")
    if (!is.numeric(exposure)) stop("exposure is not numeric.")
    if (any(exposure <= 0)) stop("Non-positive value(s) in exposure.")
    if (any(count < 0)) stop("Negative value(s) in count.")
    if (length(count) != length(exposure)) stop("length mismatch.")
    result <- cbind(count, exposure)
    colnames(result) <- c("count", "exposure")
    class(result) <- "oe"
    result
}
