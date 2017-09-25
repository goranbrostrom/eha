#' Piecewise hazards
#' 
#' Calculate piecewise hazards, no. of events, and exposure times in each
#' interval indicated by cutpoints.
#' 
#' Exact calculation.
#' 
#' @param enter Left interval endpoint
#' @param exit Right interval endpoint
#' @param event Indicator of event
#' @param cutpoints Vector of cutpoints
#' @return A list with components
#' 
#' \item{events}{Vector of number of events} \item{exposure}{Vector of total
#' exposure time} \item{intensity}{Vector of hazards, \code{intensity} ==
#' \code{events} / \code{exposure}}
#' @author Göran Broström
#' @seealso \code{\link{perstat}}
#' @keywords survival nonparametric
#' @export piecewise
piecewise <- function(enter, exit, event, cutpoints){
    n <- length(cutpoints) + 1 ## No. of time intervals.
    d <- numeric(n) ## Events
    tt <- numeric(n) ## Risk times
    
    ## assume 0 <= enter < exit < \infty.
    
    nn <- length(enter) ## Check length(exit), length(event), etc.
    
    ## First interval:
    d[1] <- sum( event[( (exit <= cutpoints[1]) & (exit > 0) )] )
    left <- pmin( enter, cutpoints[1] )
    right <- pmin( exit, cutpoints[1] )
    tt[1] <- sum(right - left)
    
    ## Intervals 2, ..., (n - 1):
    for ( j in 2:(n-1) ){
        d[j] <- sum( event[( (exit <= cutpoints[j]) &
                            (exit > cutpoints[j-1]) )] )
        left <- pmin( pmax(enter, cutpoints[j-1]), cutpoints[j])
        right <- pmax( pmin(exit, cutpoints[j]), cutpoints[j-1] )
        tt[j] <- sum(right - left)
    }
    
    ## Last interval:
    d[n] <- sum( event[ (exit > cutpoints[n - 1]) ] )
    left <- pmax( enter, cutpoints[n-1] )
    right <- pmax( exit, cutpoints[n-1] )
    tt[n] <- sum(right - left)
    
    intensity <- ifelse(tt > 0, d / tt, NA)
    list(events = d, exposure = tt, intensity = intensity)
}
