#' The Log-rank test
#' 
#' Performs the log-rank test on survival data, possibly stratified.
#' 
#' @usage logrank(formula, data = parent.frame())
#' 
#' @param formula formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function.
#' @param data a data.frame in which to interpret the variables named in the
#' formula.
#' @return A list of class \code{logrank} with components
#' \item{test.statistic}{The logrank (score) test statistic.}
#' \item{df}{The degrees of freedom of the test statistic.}
#' \item{p.value}{The p value of the test.}
#' \item{hazards}{A list of two-column matrices, describing event times and 
#' corresponding hazard atoms in each stratum (class 'hazdata').}
#' @note The test is performed by fitting a Cox regression model and reporting
#' its \code{score test}. With tied data, this might be slightly different from
#' the true logrank test, but the difference is unimprtant in practice.
#' @author Göran Broström
#' @seealso \code{\link{coxreg}}, \code{\link{print.logrank}}.
#' @keywords survival
#' @examples 
#' fit <- logrank(Surv(enter, exit, event) ~ sex, data = oldmort)
#' fit
#' @export    
logrank <- function(formula, data = parent.frame()){
    fit <- coxreg(formula, data, coxph = TRUE)
    tval <- fit$score
    df <- fit$df
    pval <-pchisq(tval, df, lower.tail = FALSE)
    fit <- list(test.statistic = tval, df = df, p.value = pval,
                hazards = fit$hazards)
    class(fit) <- "logrank"
    fit
}