#' Convert time in years since "0000-01-01" to a date.
#' 
#' This function uses \code{as.Date} and a simple linear transformation.
#' 
#' 
#' @param times a vector of durations
#' @return A vector of dates as character strings of the type "1897-05-21".
#' @author Göran Broström
#' @seealso \code{\link{toTime}}
#' @keywords survival
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' toDate(1897.357)
#' 
#' @export toDate
toDate <- function(times){
    if (!is.numeric(times)) stop("Argument must be numeric")
    times * 365.2425 + as.Date("0000-01-01")
}
