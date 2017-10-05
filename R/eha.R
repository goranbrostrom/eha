#' @details
#' 
#' Eha enhances the recommended \pkg{survival} package in several ways, 
#' see the description. The main applications in mind are demography and 
#' epidemiology. For standard Cox regression analysis the function 
#' \code{\link{coxph}} in \pkg{survival} is still recommended. The function 
#' \code{\link{coxreg}} in \pkg{eha} in fact calls coxph for the standard kind 
#' of analyses.
#' 
#' @references 
#' Broström, G. (2012). \emph{Event History Analysis with R}, Chapman and Hall/CRC 
#' Press, Boca Raton, FL.
#' @import survival
#' @import stats
#' @import graphics
#' @useDynLib eha
"_PACKAGE"
