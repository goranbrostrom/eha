#' Plots of hazdata objects.
#' 
#' Baseline hazards estimates.
#' 
#' It is also possible to have as first argument an object of type "coxreg",
#' given that it contains a component of type "hazdata".
#' 
#' @param x A \code{logrank} object, typically the 'hazards' element in the
#' output from \code{link{logrank}}.
#' @param fn Which type of plot? Allowed values are "cum" (or "cumhaz"),
#' "surv" (or "sur"), "log", or "loglog". The last two plots the cumulative 
#' hazards on a log (y) scale or a log-log (xy) scale, respectively. 
#' @param xlim Horizontal plot limits. If NULL, calculated by the function.
#' @param ylim Vertical plot limits. If NULL, set to \code{c(0, 1)} for a 
#' plot of the survival function.
#' @param main A heading for the plot.
#' @param xlab Label on the x axis.
#' @param ylab Label on the y-axis.
#' @param col Color of the lines. May be a vector of length equal to No. of
#' strata.
#' @param lty Line type(s). May be a vector of length equal to No. of strata.
#' @param printLegend Logical or character; should a legend be produced?  
#' Defaults to TRUE. If character, it should be one of \code{bottomleft}, 
#' \code{bottomright}, etc, see \code{\link{legend}}.
#' @param ... Anything that \code{\link[graphics]{plot.default}} likes...
#' @return A list where the elements are two-column matrices, one for each
#' stratum in the model. The first column contains risktimes, and the second
#' the y coordinates for the requested curve(s).
#' @note \code{x} is a list where each element is a two-column matrix. The first
#'  column contains failure times, and the second column contains the 
#'  corresponding 'hazard atoms'. 
#' @author Göran Broström
#' @keywords survival
#' @examples
#' 
#' fit <- logrank(Surv(enter, exit, event), group = civ, data = oldmort[oldmort$region == "town", ])
#' plot(fit)
#' 
#' @export
plot.logrank <- function(x,
                         fn = c("cum", "surv", "log", "loglog"),
                         xlim = NULL,
                         ylim = NULL,
                         main = NULL,
                         xlab = "",
                         ylab = "",
                         col = "black",
                         lty = 1,      # Jan 17, 2014.
                         printLegend = TRUE, # Jul 3, 2020
                         ...
                         ){
    if (!inherits(x, "logrank")) stop("Only for 'logrank' objects.")
    plot(x$hazards, fn = fn, xlim = xlim, ylim = ylim, main = main,
         xlab = xlab, ylab = ylab, col = col, lty = lty, 
         printLegend = printLegend, ...)
}

