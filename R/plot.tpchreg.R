#' Plots output from a tpchreg regression
#' 
#' Plot(s) of the hazard, cumulative hazards, and/or the survivor
#' function(s) for each stratum.
#' 
#' 
#' @param x A \code{tpchreg} object
#' @param fn Which functions shoud be plotted! Default is the hazard function.
#' @param log character, "" (default), "y", or "xy".
#' @param main Header for the plot
#' @param xlim x limits
#' @param ylim y limits
#' @param xlab x label
#' @param ylab y label
#' @param col Color(s) for the curves. Defaults to black.
#' @param lty Line type for the curve(s). Defaults to 1:(No. of strata).
#' @param printLegend Logical, or character ("topleft", "bottomleft",
#' "topright" or "bottomright"); if \code{TRUE} or character, a legend is added
#' to the plot if the number of strata is two or more.
#' @param \dots Extra parameters passed to 'plot' and 'lines'.
#' @return No return value.
#' @author Göran Broström
#' @seealso \code{\link{tpchreg}}
#' @keywords dplot survival
#' 
#' @export
plot.tpchreg <- function(x,
                       fn = c("haz", "cum", "sur"),
                       log = "",
                       main = NULL,
                       xlim = NULL,
                       ylim = NULL,
                       xlab = "Duration",
                       ylab = "",
                       col, 
                       lty, 
                       printLegend = TRUE,
                       ...){

    if (!inherits(x, "tpchreg")) stop("Works only with 'tpchreg' objects.")
    if (missing(col)) col <- rep(1, x$n.strata) ## New 2013-12-05
    if (missing(lty)) lty <- 1:x$n.strata # No. of strata

    if (length(col) < x$n.strata) col <- rep(col, x$n.strata)
    if (length(lty) < x$n.strata) lty <- rep(lty, x$n.strata)

    if (!(all(fn %in% c("haz", "cum", "sur"))))
        stop(paste(fn[1], "is an illegal value of 'fn'"))


    if (length(fn) >= 3){
        oldpar <- par(mfrow = c(2, 2))
        on.exit(par(oldpar))
    }else if (length(fn) == 2){
        oldpar <- par(mfrow = c(2, 1))
        on.exit(par(oldpar))
    }
    if (is.null(xlim)){
        xlim <- range(x$cuts)
    }
    npts <- 4999
    ns <- NROW(x$hazards)
    xx <- seq(xlim[1], xlim[2], length = npts)
    haz <- matrix(0, ncol = npts, nrow = ns)
    sur <- haz
    Haz <- haz
    ## hazard
    cuts <- x$cuts - x$cuts[1]
    cuts <- cuts[-c(1, length(cuts))]
    
    ##cuts <- cuts[-length(cuts)]
        for (i in 1:ns){
            haz[i, ] <- hpch(xx - x$cuts[1], cuts, x$hazards[i, ])
            sur[i, ] <- ppch(xx - x$cuts[1], cuts, x$hazards[i, ],
                             lower.tail = FALSE)
            Haz[i, ] <- Hpch(xx - x$cuts[1], cuts, x$hazards[i, ])
        }
        dist <- "Pcwise const"

    if ("haz" %in% fn){

        if (is.null(ylim)) {
            ylim0 <- c(0, max(haz))
        }else{
            ylim0 <- ylim
        }
        if (log != ""){
            ylim0 <- NULL # ??
            xlim <- NULL
        }
        ##if (min(p) < 1) ylim0[2] <- min(ylim0[2], max(haz[, -1]))

        if (is.null(xlab)) xlab <- "Duration"
        if (is.null(ylab)) ylab <- "Hazard"
        if (is.null(main)){
            hmain <- paste(dist, "hazard function")
        }else{
            hmain <- main
        }
        plot(xx, haz[1, ], type = "l", xlim = xlim, ylim = ylim0,
             col = col[1], lty = lty[1], log = log,
             xlab = xlab, ylab = ylab, main = hmain, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, haz[i, ], type = "l", lty = lty[i],
                      col = col[i], ...) # ', ...' added in 2.4-4
            }
        }
        ##abline(h = 0)
        ##abline(v = 0)
        if (is.character(printLegend)){
            if (!(printLegend %in% c("topleft", "bottomleft",
                                     "topright", "bottomright",
                                     "bottom", "left", "top",
                                     "right", "center"))){
                printLegend <- FALSE
                warning("Illegal value of 'printLegend'")
            }
        }
        if (is.logical(printLegend)){
            if ((ns > 1) && printLegend){
                legend(x = "topleft",  legend = x$strata, lty = lty,
                       inset = 0.001,
                       col = col)
            }
        }else{
            if ((ns > 1) && is.character(printLegend)){
                    legend(x = printLegend,  legend = x$strata,
                           lty = lty, inset = 0.001, col = col)
            }
        }
    }
    ## Cumulative hazard
    if ("cum" %in% fn){
        
        if (is.null(ylim)){
            ylim0 <- c(0, max(Haz))
        }else{
            ylim0 <- ylim
            ylim0[2] <- max(ylim0[2], max(Haz))
        }
        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Cumulative Hazard"
        if (is.null(main)){
            Hmain <- paste(dist, "cumulative hazard function")
        }else{
            Hmain <- main
        }
        plot(xx, Haz[1, ], type = "l", xlim = xlim, ylim = ylim0,
             xlab = xlab, ylab = ylab, main = Hmain, col = col[1],
             lty = lty[1], ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, Haz[i, ], type = "l", lty = lty[i],
                      col = col[i], ...) # ', ...' added in 2.4-4
            }
        }
        ##abline(h = 0)
        ##abline(v = 0)
        if (is.logical(printLegend)){
            if ((ns > 1) && printLegend){
                legend(x = "topleft",  legend = x$strata, lty = lty,
                       col = col, inset = 0.001)
            }
        }else{
            if ((ns > 1) && is.character(printLegend)){
                    legend(x = printLegend, legend = x$strata,
                           lty = lty, inset = 0.001, col = col)
            }
        }
    }

    ## Survivor function
    if ("sur" %in% fn){
        ##if (is.null(ylim))
        ylim <- c(0, 1)

        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Survival"
        if (is.null(main)){
            smain <- paste(dist, "survivor function")
        }else{
            smain <- main
        }
        plot(xx, sur[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = smain, lty = lty[1],
             col = col[1], ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, sur[i, ], type = "l", lty = lty[i],
                      col = col[i], ...) # ', ...' added in 2.4-4
            }
        }
        abline(h = 0)
        abline(v = 0)
        if (is.logical(printLegend)){
            if ((ns > 1) && printLegend){
                legend(x = "bottomleft",  legend = x$strata, lty = lty,
                       inset = 0.001,
                       col = col)
            }
        }else{
            if ((ns > 1) && is.character(printLegend)){
                    legend(x = printLegend,  legend = x$strata,
                           lty = lty, inset = 0.001, col = col)
            }
        }
    }
}

