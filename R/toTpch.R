#' Transform survival data to tabular form
#' 
#' Transform a "survival data frame" to tabular form aggregating number of events 
#' and exposure time by time intervals and covariates. 
#' 
#' @usage toTpch(formula, data, cuts, enter = "enter", exit = "exit",
#' event = "event", episode = "age")
#' 
#' @param formula A model formula.
#' @param data A data frame with survival data.
#' @param cuts An ordered, non-negative vector of time points at which a hazard function changes value.
#'  Note that data are left truncated at cuts[1] (the smallest value) and right censored at c[n], where 
#'  n  is the length of cuts and cuts[n] == max(cuts).
#' @param enter Character string with the name of the variable representing left truncation values.
#' @param exit Character string with the name of the event/censoring time variable.
#' @param event Character string with the name of the event indicator variable.
#' @param episode Character string with the name of the output variable of the grouped time (a factor variable).
#' @details If \code{cuts} is missing, nothing is done. Internally, this function first calls 
#'  \code{survival::survSplit} and the \code{stats::aggregate}. 
#' @note Episodes, or parts of episodes, outside \code{min(cuts), max(cuts)} are cut off. 
#'  With continuous covariates, consider rounding them so that the number of distinct oberved values is not too large.
#' @return A data frame with exposure time and number of events aggregated by time intervals and covariates.
#'  If all covariates are factors,this usually results in a huge reduction of the size of thedata frame, 
#'  but otherwise the size of the output may be larger than the size of the input data frame  
#' @author Göran Broström
#' @keywords manip
#'  
#' @export    
toTpch <- function(formula, data, cuts, enter = "enter", exit = "exit",
                   event = "event", episode = "age"){
    ## cuts are ordered distinct non-negative numbers
    ## 0 <= c[1] < c[2] < ... < c[n] < Infty, 
    ## defining  n - 1 intervals
    ## (c[1], c[2]], (c[2], c[3]], ..., (c[n-1], c[n]].
    ##
    ## Note that this implies that data are left truncated at c[1]
    ## (set = 0 if no left truncation), and 
    ## that the last cut, c[n], is effectively == Infinity. 
    
    if (missing(cuts)){
        stop("Argument 'cuts' is missing with no default.")
    }
    n <- length(cuts)
    if (n >= 2){
        levs <- character(n-1)
        for (i in 1:(n-1)){
            levs[i] <- paste(cuts[i], cuts[i + 1], sep = "-")
        }
    }
    form_ch <- as.character(formula)
    
    resp <- form_ch[2]
    sk <- survival::survSplit(formula, data = data, cut = cuts, 
                              episode = episode, start = enter, end = exit)
    ## Trim: 
    weq <- sk[[episode]] %in% c(1, n + 1)
    sk <- sk[!weq, ]
    ##
    
    sk$exposure <- sk$exit - sk$enter
    
    new_resp <- "cbind(event, exposure)"
    form_ch[2] <- new_resp
    form_ch[3] <- paste(form_ch[3], episode, sep = " + ")

    form <- as.formula(paste(form_ch[2], form_ch[1], form_ch[3], sep = "")) 

    ##return(form)
    out <- aggregate(form, data = sk, FUN = sum)
    ## Fungerar bra...nu måste vi få in tiden (enter, exit)
    pic <- (min(out[[episode]]) - 1):(max(out[[episode]] - 1))
    levs <- levs[pic]
    out[[episode]] <- factor(out[[episode]], labels = levs)
    out
}
