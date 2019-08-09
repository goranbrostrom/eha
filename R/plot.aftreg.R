#' Plots output from an AFT regression
#' 
#' Just a simple plot of the hazard (cumulative hazard, density, survival)
#' functions for each stratum.
#' 
#' The plot is drawn at the mean values of the covariates, by default.
#' 
#' @param x A \code{aftreg} object
#' @param fn Which functions shoud be plotted! Default is all. They will scroll
#' by, so you have to take care of explicitly what you want to be produced.
#' See, eg, \code{par(mfrow = ...)}
#' @param main Header for the plot
#' @param xlim x limits
#' @param ylim y limits
#' @param xlab x label
#' @param ylab y label
#' @param col Colors?
#' @param lty Line types?
#' @param printLegend Should legend be printed? Default is yes.
#' @param new.data At which covariate values?
#' @param \dots Extra parameters passed to 'plot'
#' @return No return value.
#' @author Göran Broström
#' @seealso \code{\link{aftreg}}
#' @keywords dplot survival
#' @examples
#' 
#' y <- rllogis(40, shape = 1, scale = 1)
#' x <- rep(c(1,1,2,2), 10)
#' fit <- aftreg(Surv(y, rep(1, 40)) ~ x, dist = "loglogistic")
#' plot(fit)
#' 
#' @export
plot.aftreg <- function(x,
                        fn = c("haz", "cum", "den", "sur"),
                        main = NULL,
                        xlim = NULL,
                        ylim = NULL,
                        xlab = "Duration",
                        ylab = "",
                        col,
                        lty,
                        printLegend = TRUE,
                        new.data = x$means,
                         ...){
    if (!inherits(x, "aftreg")) stop("Works only with 'aftreg' objects.")
    ##if (x$pfixed) stop("True exponential hazards are not plotted")
    if (!(all(fn %in% c("haz", "cum", "den", "sur"))))
        stop(paste(fn, "is an illegal value of 'fn'"))


    if (length(fn) >= 3){
        oldpar <- par(mfrow = c(2, 2))
        on.exit(par(oldpar))
    }else if (length(fn) == 2){
        oldpar <- par(mfrow = c(2, 1))
        on.exit(par(oldpar))
    }
    ncov <- length(x$means)
    ns <- x$n.strata
    param.scale <- if (x$param=="lifeAcc") -1 else 1
    if (x$pfixed){
        p <- rep(x$shape, ns)
        lambda <- exp(x$coefficients[ncov + (1:ns)] +
                      param.scale*sum((new.data - x$means) * x$coefficients[1:ncov]))
    }else{
        p <- exp(x$coefficients[ncov + (1:ns) * 2])
        lambda <- exp(x$coefficients[ncov + (1:ns) * 2 - 1] +
                      param.scale*sum((new.data - x$means) * x$coefficients[1:ncov]))
    }

    if (ncov){
        score <- exp(sum((new.data - x$means) * x$coefficients[1:ncov]))
    }else{
        score <- 1
    }

    ##if (ncov){ # THIS IS for aftplot!!
    ##    uppe <- exp(-sum((new.data[1:ncov] * x$coefficients[1:ncov]) / p)
    ##    lambda <- lambda * uppe
    ##}
    if (is.null(xlim))
        xlim <- c(min(x$y[, 1]), max(x$y[, 2]))

    npts <- 4999
    xx <- seq(xlim[1], xlim[2], length = npts)
    ##if (xx[1] <= 0) xx[1] <- 0.001

    ## hazard
    if (x$dist == "weibull"){
        dist <- "Weibull"
        haza <- hweibull
        Haza <- Hweibull
        Surviv <- pweibull
        Dens <- dweibull
    }else if (x$dist == "loglogistic"){
        dist <- "Loglogistic"
        haza <- hllogis
        Haza <- Hllogis
        Surviv <- pllogis
        Dens <- dllogis
    }else if (x$dist == "lognormal"){
        dist = "Lognormal"
        haza <- hlnorm
        Haza <- Hlnorm
        Surviv <- plnorm
        Dens <- dlnorm
    }else if (x$dist == "ev"){
        dist = "Extreme value"
        haza <- hEV
        Haza <- HEV
        Surviv <- pEV
        Dens <- dEV
    }else if (x$dist == "gompertz"){
        dist = "Gompertz"
        haza <- hgompertz
        Haza <- Hgompertz
        Surviv <- pgompertz
        Dens <- dgompertz
        ## canonical parameterisation
        p <- p/lambda
    }

    if ("haz" %in% fn){
        haz <- matrix(ncol = npts, nrow = ns)
        for (i in 1:ns){
            haz[i, ] <- haza(xx, scale = lambda[i], shape = p[i])
        }

        if (is.null(ylim)) ylim <- c(0, max(haz))
        if (min(p) < 1) ylim[2] <- min(ylim[2], max(haz[, -1]))

        if (is.null(xlab)) xlab <- "Duration"
        if (is.null(ylab)) ylab <- "Hazard"
        if (is.null(main)) main <- paste(dist, "hazard function")
        plot(xx, haz[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, ...)
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
            Haz[i, ] <- Haza(xx, scale = lambda[i], shape = p[i])
        }
        ylim <- c(0, max(Haz))
        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Cumulative Hazard"
        ##if (is.null(main))
        main <- paste(dist, "cumulative hazard function")
        plot(xx, Haz[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, ...)
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
            if (dist == "Lognormal"){
                sdlog <- 1 / p[i]
                meanlog <- log(lambda[i])
                den[i, ] <- dlnorm(xx, meanlog, sdlog)
            }else{
                den[i, ] <- Dens(xx, scale = lambda[i], shape = p[i])
            }
        }

        ##if (is.null(ylim))
        ylim <- c(0, max(den))

        if (min(p) < 1) ylim[2] <- min(max(den[, -1]))

        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Density"
        ##if (is.null(main))
        main <- paste(dist, "density function")
        plot(xx, den[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, ...)
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
            if (dist == "Lognormal"){
                sdlog <- 1 / p[i]
                meanlog <- log(lambda[i])
                sur[i, ] <- plnorm(xx, meanlog, sdlog,
                                   lower.tail = FALSE)
            }else{
                sur[i, ] <- Surviv(xx, scale = lambda[i],
                                   shape = p[i],
                                   lower.tail = FALSE)
            }
        }

        ##if (is.null(ylim))
        ylim <- c(0, 1)

        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Survival"
        ##if (is.null(main))
        main <- paste(dist, "survivor function")
        plot(xx, sur[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, ...)
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
