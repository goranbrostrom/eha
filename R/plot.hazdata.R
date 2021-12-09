#' Plots of hazdata objects.
#' 
#' Baseline hazards estimates.
#' 
#' It is also possible to have as first argument an object of type "coxreg",
#' given that it contains a component of type "hazdata".
#' 
#' @param x A \code{hazdata} object, typically the 'hazards' element in the
#' output from \code{link{coxreg}} with \code{method = "ml"} or 
#' \code{method = "mppl"} or \code{coxph = FALSE}.
#' @param strata Stratum names if there are strata present.
#' @param fn Which type of plot? Allowed values are "cum" (or "cumhaz"),
#' "surv" (or "sur"), "log", or "loglog". The last two plots the cumulative 
#' hazards on a log (y) scale or a log-log (xy) scale, respectively. 
#' @param fig Should a plot actually be produced? Default is TRUE.
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
#' time0 <- numeric(50)
#' group <- c(rep(0, 25), rep(1, 25))
#' x <- runif(50, -0.5, 0.5)
#' time1 <- rexp( 50, exp(group) )
#' event <- rep(1, 50)
#' fit <- coxreg(Surv(time0, time1, event) ~ x + strata(group), method = "ml")
#' plot(fit$hazards, col = 1:2, fn = "surv", xlab = "Duration")
#' ## Same result as:
#' ## plot(fit, col = 1:2, fn = "sur", xlab = "Duration")
#' 
#' @export
plot.hazdata <- function(x, strata = NULL,
                         fn = c("cum", "surv", "log", "loglog"),
                         fig = TRUE,
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
    ## Added 7 dec 2013: col.
    ## x is of type 'hazdata', which is a list with two-column matrices
    ## as components, one component per stratum. The first column contains
    ## risktimes, and the second column the corresponding 'hazard atoms'.
    if (!inherits(x, "hazdata")) stop("Only for 'hazdata' objects.")
    where <- "topleft"
    if (is.logical(printLegend) && printLegend){
        if (fn[1] %in% c("surv", "sur")){
            where <- "bottomleft"
        }else{
            where <- "topleft"
        }
    }else if (is.character(printLegend)){
        where <- printLegend
        printLegend <- TRUE
    }else if (!is.logical(printLegend)){
        stop("Wrong value to 'printLegend'")
    }

    if (!(where %in% c("bottomleft", "bottomright", "topleft", "topright",
                       "left", "right", "top", "bottom", "center")))
        stop(paste(where, " is not allowed as a value of 'where'"))
    
    if (is.null(x)){
        cat("Must be fixed in plot.hazdata!\n")
        return(NULL)
    }
        
    # if (!inherits(x, "hazdata")){
    #     if (!inherits(x, "coxreg")){
    #         stop("First argument must be of type 'hazdata' or 'coxreg'")
    #     }else{
    #         y <- x
    #         x <- x$hazards
    #         if (is.null(x)) stop("No 'hazards' component present")
    #     }
    # }
    fn <- fn[1]
    
    if (fn == "sur"){
        fn <- "surv"
    }else{
        if (fn == "cumhaz") fn <- "cum"
    }
    if (!(fn %in% c("cum", "surv", "log", "loglog")))
        stop(paste(fn, "is an illegal value of 'fn'"))

    n.strata <- length(x)
    if (n.strata >= 2 & is.null(strata) & !is.null(names(x))){
        strata <- names(x)
    }

    if (length(col) < n.strata) col <- rep(col, n.strata)
    if (length(lty) < n.strata) lty <- 1:n.strata
    
    yVal <- function(x){
        if (fn == "cum") return(cumsum(x))
        if (fn %in% c("log", "loglog")){
            y <- cumsum(x)
            y[y <= 0] <- NA
            return(y)
        } ##return(log(cumsum(x)))
        n <- length(x)
        s <- numeric(n)
        s[1] <- 1 - x[1]
        if (n > 1){
            for (rs in 2:n){
                s[rs] <- s[rs - 1] *(1 - x[rs])
            }
        }
        return(s)
    }

    max.x <- max(x[[1]][, 1])
    min.x <- min(x[[1]][, 1])

    max.y <- max(x[[1]][, 2]) # What else :-)
    min.y <- min(x[[1]][, 2])
    for (i in 1:n.strata){
        x[[i]][, 2] <- yVal(x[[i]][, 2])
        max.y <- max(c(max.y, x[[i]][, 2]), na.rm = TRUE)
        min.y <- min(c(min.y, x[[i]][, 2]), na.rm = TRUE)

        max.x <- max(c(max.x, x[[i]][, 1]))
        min.x <- min(c(min.x, x[[i]][, 1]))
    }

    if (is.null(xlim)) {
        xlim <- c(min.x, max.x)
        ## if (fn != "loglog") xlim[1] <- 0.0 # Why? 2.2-6
    }
    if (length(xlim) != 2) stop("'xlim' must be a vector of length two")

    if (is.null(ylim)){
        ylim <- c(min.y, max.y)
        if ( (fn == "surv") || (fn == "cum") ) ylim[1] <- 0.0
    }
    if (length(ylim) != 2) stop("'ylim' must be a vector of length two")

    if (fn != "loglog"){
        for (i in 1:n.strata){
            if (fn == "surv"){
                x[[i]] <- rbind(c(x[[i]][1, 1], 1), x[[i]])
            }else{
                x[[i]] <- rbind(c(x[[i]][1, 1], 0), x[[i]]) # Not outcommented
                                                            # 2.9.0.9200!
            }
            ##x[[i]][, 1] <- c(x[[i]][1, 1], x[[i]][, 1])
        }
    }
    if (fig){
        if (fn == "log") {
            ##ylim <- NULL
            loga <- "y"
        }else if (fn == "loglog"){
            loga <- "xy"
            ##ylim <- NULL
        }else{
            loga <- ""
        }
        
        plot(x[[1]][, 1], x[[1]][, 2], type = "s", log = loga,
             xlim = xlim, ylim = ylim, 
             col = col[1],
             xlab = xlab, ylab = ylab, main = main, lty = lty[1], ...)
        if (n.strata > 1){
            for (i in 2:n.strata){
                lines(x[[i]][, 1], x[[i]][, 2], type = "s", lty = lty[i],
                      col = col[i], ...)
            }
            if (is.null(strata)) strata <- 1:n.strata
            if (printLegend){
                legend(where, legend = strata, lty = lty,
                       col = col, inset = 0.02)
            }
        }else{
            abline(h = 0)
        }
    }
    invisible(list(x = x, fn = fn))
}
