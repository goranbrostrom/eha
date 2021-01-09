#' Graphical goodness-of-fit test
#' 
#' Comparison of the cumulative hazards functions for a semi-parametric and a
#' parametric model.
#' 
#' For the moment only a graphical comparison. The arguments \code{sp} and
#' \code{pp} may be swapped.
#' 
#' @param sp An object of type "coxreg", typically output from
#' \code{\link{coxreg}}
#' @param pp An object of type "phreg", typically output from
#' \code{\link{phreg}}
#' @param main Header for the plot. Default is distribution and "cumulative
#' hazard function"
#' @param col Line colors. should be \code{NULL} (black lines) or of length 2
#' @param lty line types.
#' @param printLegend Should a legend be printed? Default is \code{TRUE}.
#' @return No return value.
#' @author Göran Broström
#' @seealso \code{\link{coxreg}} and \code{\link{phreg}}.
#' @keywords distribution
#' @examples
#' 
#' data(mort)
#' oldpar <- par(mfrow = c(2, 2))
#' fit.cr <- coxreg(Surv(enter, exit, event) ~ ses, data = mort)
#' fit.w <- phreg(Surv(enter, exit, event) ~ ses, data = mort)
#' fit.g <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
#' dist = "gompertz")
#' fit.ev <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
#' dist = "ev")
#' check.dist(fit.cr, fit.w)
#' check.dist(fit.cr, fit.g)
#' check.dist(fit.cr, fit.ev)
#' par(oldpar)
#' 
#' @export
check.dist <- function(sp, pp, main = NULL, col = 1:2, 
                       lty = 1:2, printLegend = TRUE){
    if (!inherits(sp, "coxreg")){
        if (inherits(pp, "coxreg")){ # swap:
            tmp <- pp
            pp <- sp
            sp <- tmp
            rm(tmp)
        }else{
            stop ("Some argument must be of type 'coxreg'")
        }
    }
    if (!inherits(pp, c("phreg", "pchreg", "tpchreg")))
        stop ("Some argument must be of type 'phreg' or 'pchreg' or 'tpchreg'.")

    if (length(col) == 1) col <- c(col, col)
    ##if (!sp$nullModel){ # NOTE: 'center' is deprecated in both!
      ##  if ((!sp$center) && pp$center)
        ##    warning("The non-parametric fit is not centered.") 
    ##    if ((!pp$center) && sp$center)
      ##      warning("The parametric fit is not centered.")
    ##}
    if ((!is.null(sp$strata)) || (!is.null(pp$strata)))
        stop("Not for stratified fits; try a comparison stratum by stratum.") 
    if (is.null(main)){
        main <- pp$dist # Capitalize:
        substr(main, 1, 1) <- toupper(substr(main, 1, 1))
        if (main == "Pch") main <- "Piecewise constant"
        if (main == "Ev") main = "Extreme value"
    }

    ##x.max <- max(pp$y[, 2])
    plot.coxreg(sp, fn = "cum", main = main, col = col[1], lty = lty[1])
    
    ##x <- plot(pp, fn = "cum", fig = FALSE) # Gives error with gompertz
    oo <- hazards(pp)
    lines(oo$x, oo$y, col = col[2], lty = lty[2])
    
    if (printLegend){
        legd <- c("Non-parametric", pp$dist)    
        legend(x = "topleft", legend = legd,
               col = col, lty = lty)
    }
}

