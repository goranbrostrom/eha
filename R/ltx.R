#' LaTeX printing of coxreg results.
#' 
#' This (generic) function prints the LaTeX code of the results of a fit from
#' \code{\link{coxreg}}, \code{\link{phreg}}, or \code{\link{aftreg}}, similar
#' to what \code{xtable} does for fits from other functions.
#' 
#' The result is a printout which is (much) nicer than the standard printed
#' output from \code{glm} and friends,
#' 
#' @param x The output from a call to \code{coxreg}, \code{coxreg}, or
#' \code{aftreg}
#' @param caption A suitable caption for the table.
#' @param label A label used in the LaTeX code.
#' @param dr Output from a \code{drop1} call.
#' @param digits Number of digits to be printed.
#' @param \dots Not used.
#' @return LaTeX code version of the results from a run with
#' \code{\link{coxreg}}, \code{\link{phreg}}, or \code{\link{aftreg}}.
#' @note There is no method in \code{xtable} for \code{coxreg}.
#' @author Göran Broström.
#' @seealso \code{xtable}, \code{\link{coxreg}}
#' @keywords printing
#' @examples
#' 
#' data(oldmort)
#' fit <- coxreg(Surv(enter, exit, event) ~ civ + sex, data = oldmort)
#' dr <- drop1(fit, test = "Chisq")
#' ltx(fit, dr = dr, caption = "A test example.", label = "tab:test1") 
#' 
#' @export ltx
ltx <- function(x,
                caption = NULL,
                label = NULL,
                dr = NULL,
                digits = max(options()$digits - 4, 3), ...){
    UseMethod("ltx")
}
