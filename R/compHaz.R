#' Graphical comparison of cumulative hazards
#' 
#' Comparison of the estimated baseline cumulative hazards functions for two survival models.
#' 
#' @param fit1 An object of type "coxreg", "phreg", or other output from
#' from survival fitters.
#' @param fit2 An object of type "coxreg", "phreg", or other output from
#' survival fitters.
#' @param main Header for the plot. Default is \code{NULL}.
#' @param lty line types.
#' @param col Line colors. should be \code{NULL} (black lines) or of length 2.
#' @param printLegend Should a legend be printed? Default is \code{TRUE}.
#' @return No return value.
#' @author Göran Broström
#' @seealso \code{\link{hazards}}, \code{\link{coxreg}}, and \code{\link{phreg}}.
#' @keywords distribution
#' @examples
#' 
#' fit.cr <- coxreg(Surv(enter, exit, event) ~ sex, data = oldmort)
#' fit.w <- phreg(Surv(enter, exit, event) ~ sex, data = oldmort)
#' compHaz(fit.cr, fit.w)
#' @export
compHaz <- function(fit1, fit2, main = NULL, lty = 1:2, col = c("red", "blue"), printLegend = TRUE){
    if (fit1$n.strata > 1 || fit2$n.strata > 1) stop("Not for stratified fits.")
    h1 <- hazards(fit1)
    h2 <- hazards(fit2)
    op <- par(las = 1)
    on.exit(par(op))
    if (length(h1) == 1){
        x1 <- h1[[1]][, 1]
        y1 <- h1[[1]][, 2]
        fit1$dist <- "Non-parametric"
    }else{
        x1 <- h1$x
        y1 <- h1$y
    }
####
    if (length(h2) == 1){
        x2 <- h2[[1]][, 1]
        y2 <- h2[[1]][, 2]
        fit2$dist <- "Non-parametric"
    }else{
        x2 <- h2$x
        y2 <- h2$y
    }
    plot(x1, y1, type = "l", col = col[1], lty = lty[1], ylab = "Cum. hazards", 
         xlab = "Time", main = main, ylim = c(0, max(y1, y2)), xlim = c(min(x1, x2), max(x1, x2)))    
    lines(x2, y2, col = col[2], lty = lty[2])
    abline(h = 0)
    if (printLegend){
        lgd <- c(fit1$dist, fit2$dist)
       legend("topleft", legend = lgd, lty = lty, col = col)
    }
}