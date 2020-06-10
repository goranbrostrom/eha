#' @export
ltx.tpchreg <-function(x, caption = NULL, label = NULL, dr = NULL, 
                       digits=max(options()$digits - 4, 3), ...){
    ltx.coxreg(x, caption, label, dr, digits, ...)
}
