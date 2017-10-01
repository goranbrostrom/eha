#' Plot method for \code{coxreg} objects
#' 
#' A plot of a baseline function of a \code{coxreg} fit is produced, one curve
#' for each stratum.
#' 
#' @param x A \code{coxreg} object
#' @param fn What should be plotted? Default is "cumhaz", and the other choices
#' are "surv", "log", and "loglog".
#' @param fig logical. If \code{TRUE} the plot is actually drawn, otherwise
#' only the coordinates of the curve(s) are returned.
#' @param xlim Start and end of the x axis.
#' @param ylim Start and end of the y axis.
#' @param main A headline for the plot
#' @param xlab Label on the x axis.
#' @param ylab Label on the y axis.
#' @param col Color of the curves. Defaults to 'black'.
#' @param lty Line type(s).
#' @param printLegend Either a logical or a text string; if \code{TRUE}, a
#' legend is printed at a default place, if \code{FALSE}, no legend is printed.
#' Otherwise, if a text string, it should be one of "bottomleft",
#' "bottomright", "topleft", etc., see \code{\link{legend}} for all possibe
#' choices.
#' @param newdata Not used
#' @param ... Other parameters to pass to the plot.
#' @return An object of class \code{hazdata} containing the coordinates of the
#' curve(s).
#' @export
plot.coxreg <- function(x,
                        fn = c("cum", "surv", "log", "loglog"),
                        fig = TRUE,
                        xlim = NULL,
                        ylim = NULL,
                        main = NULL,
                        xlab = "Duration",
                        ylab = "",
                        col,
                        lty, 
                        printLegend = TRUE,
                        newdata = NULL,
                        ...){
   if (!inherits(x, c("coxreg"))) stop("Works only with 'coxreg' objects")
   if (!is.null(x$hazards)){
       y <- x$hazards
   }else{
       if (!is.null(x$stratum)){
           if(x$nullModel){
               y <- with(x, getHaz(y, stratum, rep(1, NROW(y))))
           }else{
               y <- with(x, getHaz(y, stratum, exp(linear.predictors)))
           }
       }else{
           if (x$nullModel){
               y <- with(x, getHaz(y, rep(1, NROW(y)), rep(1, NROW(y))))
           }else{
               y <- getHaz(x$y, rep(1, length(x$linear.predictors)),
                           exp(x$linear.predictors))
           }
       }
   }
    if (missing(col)) col <- "black"
   if (missing(lty)) lty <- 1
   if (fig){
       if (is.logical(printLegend)){
           where <- NULL
       }else{
           where <- printLegend
           printLegend <- TRUE
       }
       plot.hazdata(y, strata = x$strata, fn = fn, fig = fig,
                    xlim = xlim, ylim = ylim, main = main,
                    xlab = xlab, ylab = ylab, lty = lty, col = col,
                    printLegend = printLegend, where = where, ...)
   }
   invisible(y)
}
