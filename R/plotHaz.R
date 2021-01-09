#' Graphical comparing of cumulative hazards
#' 
#' Comparison of the cumulative hazards functions for a semi-parametric and
#' parametric models.
#' 
#' For the moment only a graphical comparison. The arguments \code{sp} and
#' \code{pp} may be swapped.
#' 
#' @param sp An object of type "coxreg" or "phreg", typically output from
#' \code{\link{coxreg}} or \code{\link{phreg}}.
#' @param pp An object of type "coxreg" or "phreg", typically output from
#' \code{\link{coxreg}} or \code{\link{phreg}}.
#' @param interval Time interval for the plot, if missing, calculated from \code{sp}.
#' @param main Header for the plot. Default is distribution and "cumulative
#' hazard function"
#' @param xlab Label on x axis (default "Time")
#' @param ylab Label on y axis (default "Cum. Hazards")
#' @param col Line colors. should be \code{NULL} (black lines) or of length 2
#' @param lty line types.
#' @param ylim Y limits for the plot.
#' @param printLegend Should a legend be printed? Default is \code{TRUE}.
#' @return No return value.
#' @author Göran Broström
#' @seealso \code{\link{check.dist}}, \code{\link{coxreg}} and \code{\link{phreg}}.
#' @keywords distribution
#' @examples
#' 
#' data(mort)
#' op <- par(mfrow = c(1, 2))
#' fit.cr <- coxreg(Surv(enter, exit, event) ~ ses, data = mort)
#' fit.w <- phreg(Surv(enter, exit, event) ~ ses, data = mort)
#' fit.g <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
#' dist = "gompertz")
#' plotHaz(fit.cr, fit.w, interval = c(0, 20), main = "Weibull")
#' plotHaz(fit.cr, fit.g, main = "Gompertz")
#' par(op)
#' 
#' @export
plotHaz <- function(sp, pp, interval, main = NULL, xlab = "Time", ylab = "Cum. hazards",
                    col = c("blue", "red"), lty = 1:2, ylim, printLegend = TRUE){
   fpcl <- c("phreg", "pchreg", "tpchreg")
   if (!(inherits(sp, "coxreg") | inherits(sp, fpcl))){
      stop("Wrong first argument")
   }
   if (!(inherits(pp, "coxreg") | inherits(pp, fpcl))){
      stop("Wrong second argument")
   }
   
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
      main <- ""
   }
   
   if (missing(interval)){

      Hsp <- hazards(sp)
      Hpp <- hazards(pp)
   }else{
      Hsp <- hazards(sp, ivl = interval)
      Hpp <- hazards(pp, ivl = interval)
   }
   if (inherits(Hsp, "parhazdata")){
      x1 <- Hsp$x
      y1 <- Hsp$y
      name1 <- sp$dist
   }else{
      x1 <- Hsp[[1]][, 1]
      y1 <- Hsp[[1]][, 2]
      name1 <- "nelson-aalen"
   }
   
   if (inherits(Hpp, "parhazdata")){
      x2 <- Hpp$x
      y2 <- Hpp$y
      name2 <- pp$dist
   }else{
      x2 <- Hpp[[1]][, 1]
      y2 <- Hpp[[1]][, 2]
      name2 <- "nelson-aalen"
   }
   
   if (missing(interval)){
      interval <- c(min(x1, x2), max(x1, x2))   
   }
   
   if (missing(ylim)){
       yint <- c(0, max(y1, y2))
   }else{
      yint <- ylim
   }
   plot(x1, y1, type = "l", col = col[1], lty = lty[1], xlim = interval, ylim = yint, 
        xlab = xlab, ylab = ylab, lwd = 1.5)
   lines(x2, y2, col = col[2], lty = lty[2], lwd = 1.5)
   abline(h = 0)
   if (printLegend){
      legend(x = "topleft", legend = c(name1, name2),
             col = col, lty = lty)
   }
}

