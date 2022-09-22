#' Prints a summary of the content of a set of risk sets.
#' 
#' Given the output from \code{risksets}, summary statistics are given for it.
#' 
#' 
#' @param x An object of class 'risksets'.
#' @param ... Not used for the moment.
#' @return No value is returned; the function prints summary statistics of risk
#' sets. 
#' @note There is no \code{summary.risksets} yet. On the TODO list.
#' @author Göran Broström
#' @seealso \code{risksets}
#' @keywords summary
#' @examples
#' 
#' rs <- with(mort, risksets(Surv(enter, exit, event)))
#' print(rs)
#' 
#' @export
print.risksets <- function(x, ...){
##    if (class(x) != "risksets") stop("Only for class 'risksets'")
    if (!inherits(x, "risksets")) stop("Only for class 'risksets'")
    cat("No of strata: ", x$ns, "\n")
    if (x$ns > 1){
        for (i in 1:x$ns){
            cat("Stratum No. ", i, "\n")
            cat("--------------\n")
            cat("Number of risksets: ", x$antrs[i], "\n")
            cat("Number of events per risk set: ", summary(x$n.events), "\n")
        }
    }
}
