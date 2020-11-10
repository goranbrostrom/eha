#' The Log-rank test
#' 
#' Performs the log-rank test on survival data, possibly stratified.
#' 
#' @usage logrank(Y, group, data = parent.frame())
#' 
#' @param Y a survival object as returned by the \code{\link{Surv}} function.
#' @param group defines the groups to be compared. Coerced to a factor.
#' @param data a data.frame in which to interpret the variables.
#' @return A list of class \code{logrank} with components
#' \item{test.statistic}{The logrank (score) test statistic.}
#' \item{df}{The degrees of freedom of the test statistic.}
#' \item{p.value}{The p value of the test.}
#' \item{hazards}{A list of two-column matrices, describing event times and 
#' corresponding hazard atoms in each stratum (class 'hazdata').}
#' \item{call}{The call}
#' @note The test is performed by fitting a Cox regression model and reporting
#' its \code{score test}. With tied data, this might be slightly different from
#' the true logrank test, but the difference is unimportant in practice.
#' @author Göran Broström
#' @seealso \code{\link{coxreg}}, \code{\link{print.logrank}}.
#' @keywords survival
#' @examples 
#' fit <- logrank(Y = Surv(enter, exit, event), group = civ, data = oldmort)
#' fit
#' @export    
logrank <- function(Y, group, data = parent.frame()){
    cl <- match.call()
    Y1 <- eval(substitute(Y), data)
    g1 <- eval(substitute(group), data)
    X <- model.matrix(~ g1)[, -1, drop = FALSE]
    fit <- coxreg.fit(X, Y1, max.survs = NROW(Y1))

    tval <- fit$score
    df <- fit$df
    pval <-pchisq(tval, df, lower.tail = FALSE)
    hazards <- getHaz(Y1, strats = g1)
    fit <- list(test.statistic = tval, df = df, p.value = pval,
                hazards = hazards, call = cl)
    class(fit) <- "logrank"
    fit
}