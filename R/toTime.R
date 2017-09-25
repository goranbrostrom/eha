#' Calculate duration in years from "0000-01-01" to a given date
#' 
#' Given a vector of dates, the output is a vector of durations in years since
#' "0000-01-01".
#' 
#' 
#' @param dates A vector of dates in character form or of class \code{Date}
#' @return A vector of durations, as decribed above.
#' @author Göran Broström
#' @seealso \code{\link{toDate}}
#' @keywords survival
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' toTime(c("1897-05-16", "1901-11-21"))
#' 
#' @export toTime
toTime <- function(dates){
    if (is.numeric(dates)) return(dates)
    dates <- as.character(dates)
    c(as.Date(dates) - as.Date("0000-01-01") ) / 365.2425
  }
