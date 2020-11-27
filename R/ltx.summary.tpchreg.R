#' @export
ltx.summary.tpchreg <-function(x, caption = NULL, label = NULL, dr = NULL, 
                       digits=max(options()$digits - 4, 3), ...){
    if (!inherits(x, "summary.tpchreg")){
        stop("Wrong object type.")
    }
    class(x) <- c(class(x), "coxreg")
    ltx.coxreg(x, caption, label, dr = x$dr, digits, ...)
}
