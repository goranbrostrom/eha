#' Get baseline hazards atoms from fits from 
#' 
#' @param x A \code{reg} object.
#' @param cum Logical: Should the cumulative hazards be returned?
#' @param ... Additional arguments for various methods.
#' @return A list where each component is a two-column matrix representing
#' hazard atoms from one stratum. The first column is event time, and the 
#' second column is the corresponding hazard atom.
#' @export
hazards <- function(x, cum = TRUE, ...){
    ## x is output from '...reg'.
    UseMethod("hazards")
}
#' @export
hazards.coxreg <- function(x, cum = TRUE, ...){
    ## x is output from 'coxreg'.
    if (!is.null(x$hazards)){
        names(x$hazards) <- x$strata
         return(x$hazards)
    }
    if (is.null(x$y)) stop("Argument is lacking an y component.")
    ##
    haz <- getHaz(x$y, strats = x$stratum, score = exp(x$linear.predictors))
    ns <- length(haz)
    if (cum){
        for (i in 1:ns){
            haz[[i]][, 2] <- cumsum(haz[[i]][, 2])
        }
    }
    if (ns > 1){
        names(haz) <- x$strata
    }
    haz
}
#' @export
hazards.phreg <- function(x, cum = TRUE, ivl, n.points = 300, ...){
    if (missing(ivl)){
        if (is.null(x$y))stop("No time vector and no y value, either needed.")
        ivl <- numeric(2)
        if (ncol(x$y) == 2){
            ivl[1] <- 0
            ivl[2] <- max(x$y[, 1])
        }else{
            ivl[1] <- min(x$y[, 1])
            ivl[2] <- max(x$y[, 2])
        }
    }
    xp <- seq(ivl[1], ivl[2], length = n.points)
    if (x$dist == "pch"){
        ns <- nrow(x$hazards)
        yp <- matrix(0, nrow = ns, ncol = n.points)
        for (i in 1:ns){
            if (cum){
                yp[i, ] <- Hpch(xp, cuts = x$cuts, levels = x$hazards[i, ])
            }else{
                yp[i, ] <- hpch(xp, cuts = x$cuts, levels = x$hazards[i, ])
            }
        }
    }else{
        ns <- x$n.strata
        yp <- matrix(0, nrow = ns, ncol = n.points)
        scsh <- x$coefficients[-(1:x$df)] # Has length 2 * ns
        sc <- numeric(ns)
        sh <- numeric(ns)
        for (i in 1:ns){
            sc[i] <- scsh[2 * i - 1]
            
            sh[i] <- exp(scsh[2 * i])
            ##cat("i = ", sc[i], sh[i], "\n")
        }
        ##if (x$dist != "gompertz" |
          ##  (x$dist == "gompertz" & x$param != "rate")){
            ##sc <- exp(sc)
        ##}
        if (x$dist == "weibull"){
            sc <- exp(sc)
            for (i in 1:ns){
                if (cum){
                    yp[i, ] <- Hweibull(xp, scale = sc[i], shape = sh[i])
                }else{
                    yp[i, ] <- hweibull(xp, scale = sc[i], shape = sh[i])
                }
                
            }
        }
        if (x$dist == "ev"){
            sc <- exp(sc)
            for (i in 1:ns){
                if (cum){
                    yp[i, ] <- HEV(xp, scale = sc[i], shape = sh[i])
                }else{
                    yp[i, ] <- hEV(xp, scale = sc[i], shape = sh[i])
                }
                
            }
        }
        if (x$dist == "gompertz"){
            if (x$param == "canonical"){
                sc <- exp(sc)
                for (i in 1:ns){
                    if (cum){
                        yp[i, ] <- Hgompertz(xp, scale = sc[i], shape = sh[i], 
                                             param = "canonical")
                    }else{
                        yp[i, ] <- hgompertz(xp, scale = sc[i], shape = sh[i],
                                             param = "canonical")
                    }
                }  
            }
            if (x$param == "rate"){
                for (i in 1:ns){
                    if (cum){
                        yp[i, ] <- Hgompertz(xp, rate = sc[i], shape = sh[i], 
                                             param = "rate")
                    }else{
                        yp[i, ] <- hgompertz(xp, rate = sc[i], shape = sh[i],
                                             param = "rate")
                    }
                    
                }
            }
        }
    }
    if (!is.null(x$strata)){
        rownames(yp) <- x$strata
    }
    ret <- list(x = xp, y = yp)
    class(ret) <- "parhazdata"
    ret
}