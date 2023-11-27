#' Plots output from a Weibull regression
#' 
#' Plot(s) of the hazard, density, cumulative hazards, and/or the survivor
#' function(s) for each stratum.
#' 
#' The plot is drawn at the mean values of the covariates.
#' 
#' @param x A \code{weibreg} object
#' @param fn Which functions shoud be plotted! Default is all. They will scroll
#' by, so you have to take care explicitely what you want to be produced. See,
#' eg, \code{par(mfrow = ...)}
#' @param main Header for the plot
#' @param xlim x limits
#' @param ylim y limits
#' @param xlab x label
#' @param ylab y label
#' @param new.data At which covariate values?
#' @param \dots Extra parameters passed to 'plot'
#' @return No return value
#' @author Göran Broström
#' @seealso \code{\link{phreg}}, \code{\link{weibreg}}
#' @keywords dplot survival
#' @examples
#' 
#' y <- rweibull(4, shape = 1, scale = 1)
#' x <- c(1,1,2,2)
#' fit <- weibreg(Surv(y, c(1,1,1,1)) ~ x)
#' plot(fit)
#' 
#' @export
plot.weibreg <- function(x,
                         fn = c("haz", "cum", "den", "sur"),
                         main = NULL,
                         xlim = NULL,
                         ylim = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         new.data = x$means,
                         ...){
    if (!inherits(x, "weibreg")) stop("Works only with 'weibreg' objects.")
    ##if (x$pfixed) stop("True exponential hazards are not plotted")
    if (!(all(fn %in% c("haz", "cum", "den", "sur"))))
        stop(paste(fn, "is an illegal value of 'fn'"))

    if (length(fn) == 4){
        oldpar <- par(mfrow = c(2, 2))
        on.exit(par(oldpar))
    }
    ncov <- length(x$means)
    ns <- x$n.strata
    if (x$pfixed){
        p <- rep(x$shape, ns)
        lambda <- exp(x$coefficients[ncov + (1:ns)])
    }else{
        p <- exp(x$coefficients[ncov + (1:ns) * 2])
        lambda <- exp(x$coefficients[ncov + (1:ns) * 2 - 1])
    }
    if (ncov){
        uppe <- exp(-sum(new.data[1:ncov] * x$coefficients[1:ncov]) / p)
        lambda <- lambda * uppe
    }
    if (is.null(xlim)) {
        xlim0 <- c(min(x$y[, 1]), max(x$y[, 2]))
    }else{
        xlim0 <- xlim
    }
    
    npts <- 999
    xx <- seq(xlim0[1], xlim0[2], length = npts)
    ##if (xx[1] <= 0) xx[1] <- 0.001


    ## hazard
    if ("haz" %in% fn){
        haz <- matrix(ncol = npts, nrow = ns)
        for (i in 1:ns){
            haz[i, ] <- hweibull(xx, scale = lambda[i], shape = p[i])
        }

        if (is.null(ylim)){
            ylim0 <- c(0, max(haz))
        }else{
            ylim0 <- ylim
        }
        ## xlim0 is set globally above 
        if (min(p) < 1) ylim[2] <- min(ylim[2], max(haz[, -1]))

        if (is.null(xlab)){
             xlab0 <- "Duration"
        }else{
            xlab0 <- xlab
        }
        if (is.null(ylab)) {
            ylab0 <- "Hazard"
        }else{
            ylab0 <- ylab
        }
        main0 <- ifelse(is.null(main), "Weibull hazard function", main)
        plot(xx, haz[1, ], type = "l", xlim = xlim0, ylim = ylim0,
             xlab = xlab0, ylab = ylab0, main = main0, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, haz[i, ], type = "l", lty = i)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    ## Cumulative hazard
    if ("cum" %in% fn){

        Haz <- matrix(ncol = npts, nrow = ns)

    ##if (is.null(ylim))
        for (i in 1:ns){
            Haz[i, ] <- Hweibull(xx, scale = lambda[i], shape = p[i])
        }
        if(is.null(ylim)){
            ylim0 <- c(0, max(Haz))
        }else{
            ylim0 <- ylim
        }
        ##if (is.null(xlab))
        xlab0 <- ifelse(is.null(xlab), "Duration", xlab)
        ##if (is.null(ylab))
        ylab0 <- ifelse(is.null(ylab), "Cumulative Hazard", ylab)
        ##if (is.null(main))
        main0 <- ifelse(is.null(main), "Weibull cumulative hazard function", main)
        plot(xx, Haz[1, ], type = "l", xlim = xlim0, ylim = ylim0,
             xlab = xlab0, ylab = ylab0, main = main0, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, Haz[i, ], type = "l", lty = i)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    ## density
    if ("den" %in% fn){

        den <- matrix(ncol = npts, nrow = ns)
        for (i in 1:ns){
            den[i, ] <- dweibull(xx, scale = lambda[i], shape = p[i])
        }

        ##if (is.null(ylim))
        if (is.null(ylim)){
            ylim0 <- c(0, max(den))
        }else{
            ylim0 <- ylim
        }

        ##if (min(p) < 1) ylim0[2] <- min(max(den[, -1])) # What?

        ##if (is.null(xlab))
        xlab0 <- ifelse(is.null(xlab), "Duration", xlab)
        ##if (is.null(ylab))
        ylab0 <- ifelse(is.null(ylab), "Density", ylab)
        ##if (is.null(main))
        main0 <- ifelse(is.null(main), "Weibull density function", main)
        plot(xx, den[1, ], type = "l", xlim = xlim0, ylim = ylim0,
             xlab = xlab0, ylab = ylab0, main = main0, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, den[i, ], type = "l", lty = i)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    ## Survivor function
    if ("sur" %in% fn){
        sur <- matrix(ncol = npts, nrow = ns)
        for (i in 1:ns){
            sur[i, ] <- pweibull(xx, scale = lambda[i], shape = p[i],
                                 lower.tail = FALSE)
        }

        ##if (is.null(ylim))
        ylim0 <- if (is.null(ylim)) c(0, 1) else ylim

        ## if (is.null(xlab))
        xlab0 <- if (is.null(xlab)) "Duration" else xlab
        ##if (is.null(ylab))
        ylab0 <- if(is.null(ylab)) "Survival" else ylab
        ##if (is.null(main))
        main0 <- if(is.null(main)) "Weibull survivor function" else main
        plot(xx, sur[1, ], type = "l", xlim = xlim0, ylim = ylim0,
             xlab = xlab0, ylab = ylab0, main = main0, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, sur[i, ], type = "l", lty = i)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    ##par(oldpar)
}
