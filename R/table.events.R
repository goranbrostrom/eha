#' Calculating failure times, risk set sizes and No. of events in each risk set
#' 
#' From input data of the 'interval' type, with an event indicator, summary
#' statistics for each risk set (at an event time point) are calculated.
#' 
#' 
#' @param enter Left truncation time point.
#' @param exit End time point, an event or a right censoring.
#' @param event Event indicator.
#' @param strict If TRUE, then tabulating is not done after a time point where
#' all individuals in a riskset failed.
#' @return A list with components \item{times}{Ordered distinct event time
#' points.} \item{events}{Number of events at each event time point.}
#' \item{riskset.sizes}{Number at risk at each event time point.}
#' @author Göran Broström
#' @seealso \code{\link{risksets}}
#' @keywords survival
#' @examples
#' 
#' exit = c(1,2,3,4,5)
#' event = c(1,1,0,1,1)
#' table.events(exit = exit, event = event)
#' 
#' @export table.events
table.events <- function(enter = rep(0, length(exit)),
                         exit,
                         event,
                         strict = TRUE)
{
  n <- length(exit)

  ## Check input data:
  if ( length(enter) != n ) stop("enter and exit must have equal length.")
  if ( length(event) != n ) stop("event and exit must have equal length.")
  ##
  
  event <- (event != 0) ## 0 (FALSE) = censoring, else (TRUE) = event

  times <- c(unique(sort(exit[event])))
  nn <- length(times)

  rs.size <- double(nn)
  n.events <- double(nn)

  for (i in 1:nn) ## Try to avoid this loop!
    {
      rs.size[i] <- sum( (enter < times[i]) &
                        (exit >= times[i]) )
      n.events[i] <- sum( (exit == times[i]) & event )
    }

  stop.at <- which(rs.size == n.events)
  if (strict & length(stop.at))
    {
      stop.at <- min(stop.at) - 1
      if (stop.at <= 0)
          stop("First risk set is all events! Try 'strict = FALSE'")
      times <- times[1:stop.at]
      n.events <- n.events[1:stop.at]
      rs.size <- rs.size[1:stop.at]
    }
      
  return ( list(times         = times,
                events        = n.events,
                riskset.sizes = rs.size)
          )
}

