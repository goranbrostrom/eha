#' Proportional hazards regression with piecewise constant hazards and tabular 
#' data.
#' 
#' @usage tpchreg(formula, data, time, weights, last, subset, na.action, 
#' contrasts = NULL, start.coef = NULL, 
#' control = list(epsilon = 1.e-8, maxit = 200, trace = FALSE))
#' 
#' @param formula a formula with 'oe(count, exposure) ~ x1 + ...'
#' @param data a data frame with occurrence/exposure data plus covariates.
#' @param time the time variable, a factor character vector indicating time 
#' intervals, or numeric, indicating the start of intervals.
#' @param weights Case weights.
#' @param last If \code{time} is numeric, the closing of the last interval.
#' @param subset subset of data, not implemented yet.
#' @param na.action Not implemented yet.
#' @param contrasts Not implemented yet.
#' @param start.coef For the moment equal to zero, not used.
#' @param control list of control parameters for the optimization. 

#' @note The interpretation of cuts is different from that in \code{\link{hpch}}.
#' This is intentional.
#' 
#' @seealso \code{\link{oe}}.
#' @keywords survival regression table
#' @examples 
#' 
#' sw <- swepop
#' sw$deaths <- swedeaths$deaths
#' fit <- tpchreg(oe(deaths, pop) ~ strata(sex) + I(year - 2000), 
#' time = age, last = 101, data = sw[sw$year >= 2000, ])
#' summary(fit)
#' 
#' @export
tpchreg <- function(formula,
                    data,
                    time,
                    weights,
                    last,
                    subset,
                    na.action,
                    contrasts = NULL,
                    start.coef = NULL,
                    control = list(epsilon = 1.e-8,
                                   maxit = 200, trace = FALSE)){
    if (missing(time)){
        message("Argument 'time' is missing, indicating a constant baseline hazard.")
        t_miss <- TRUE
    }else{
        t_miss <- FALSE
    }
    

    if (is.list(control)) {
        if (is.null(control$epsilon))
            control$epsilon <- 1e-08
        if (is.null(control$maxit))
            control$maxit <- 200
        if (is.null(control$trace))
            control$trace <- FALSE
    }
    else {
        stop("control must be a list")
    }
    #
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "time", "offset", "weights", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    special <- "strata"
    Terms <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    if (!inherits(Y, "oe")) 
        stop("Response must be an oe object")
### Y:
    count <- Y[, 1]
    exposure <- Y[, 2]
### is the response
    attr(Terms, "intercept") <- 1
    strats <- attr(Terms, "specials")$strata
    dropx <- NULL
    if (length(strats)) {
        temp <- survival::untangle.specials(Terms, "strata", 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars) == 1) 
            strata.keep <- m[[temp$vars]]
        else strata.keep <- survival::strata(m[, temp$vars], shortlabel = TRUE)
        strats <- as.numeric(strata.keep)
    }
    if (length(dropx)) 
        newTerms <- Terms[-dropx]
    else newTerms <- Terms
    X <- model.matrix(newTerms, m)
    assign <- lapply(survival::attrassign(X, newTerms)[-1], 
                     function(x) x - 1)
    X <- X[, -1, drop = FALSE]
    
    if (length(dropx)) {
        covars <- names(m)[-c(1, (dropx + 1))]
    }
    else {
        covars <- names(m)[-1]
    }
    ## Added 5 jun 2020:
    jj <- which(covars == "(time)")
    covars <- covars[-jj]
    ## End Add
    offset <-  attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0)
        rep(0, length(count))
    else if (tt == 1)
        m[[offset]]
    else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
    
    if (t_miss){
        time <- rep("0-1", length(count))
    }else{
        time <- m$"(time)"
    }
    weights <- m$"(weights)"
    if (is.null(weights)){
        weights <- rep(1, length(count))
    }
    if(is.null(strats)){
        stratum <- rep(1, length(count))
        strata <- NULL
    }else{
        stratum <- strats
        strata <- levels(as.factor(strata.keep))
    }

    ##if (is.null(time)) time <- rep("0-1", length(count))
    if (is.character(time)) {
        time <- as.factor(time)
    }else if (is.numeric(time)){ # Start points of intervals.
        cuts <- sort(unique(time))
        n <- length(cuts)
        if (missing(last)){ 
            warning("'last' is missing, needed when 'time' is numeric. Created.")
            last <- cuts[n] + 1
        }
        cuts <- c(cuts, last)
        tlev <- paste(cuts[-(n+1)], cuts[-1], sep = "-")
        time <- factor(time, levels = cuts[-(n+1)], labels = tlev)
    }

    
    if (!is.factor(time)){
        stop("Argument 'time' must be a character or factor variable")
    }
    if (any(is.na(time))){
        cat("Bad values in 'time', returned\n")
        return(time)
    }
    cuts <- as.numeric(unique(unlist(strsplit(levels(time), "-"))))
    if (any(is.na(cuts))) stop("Bad values in the 'time' variable.")
    ##n.ivls <- length(unique(time)) # not used?
   ## Calling the work horse:
    fit <- tpchreg.fit(X, count, exposure, offset, weights, stratum, time, 
                       control)
    ##Terms
    fit$covars <- covars
    fit$terms <- Terms
    ##fit$newTerms <- newTerms
    stru <- struct(m, covars, dropx, X, exposure) 
    ## 'struct' is a function in 'eha'!
    fit$isI <- stru$isI
    fit$isF <- stru$isF
    fit$isO <- stru$isO
    fit$levels <- stru$levels
    fit$call <- call
    fit$w.means <- stru$w.means
    fit$n.events <- sum(count)
    fit$ttr <- stru$ttr
    fit$n <- length(count)
    fit$cuts <- cuts
    fit$strata <- strata
    if (fit$n.strata > 1){
        if (control$trace){
            cat("levels(time) = ", levels(time), "\n")
            cat("strata = ", fit$strata, "\n")
            cat("dim(hazards) = ", dim(fit$hazards), "\n")
        }
        colnames(fit$hazards) <- levels(time)
        rownames(fit$hazards) <- fit$strata
    }else{
        names(fit$hazards) <- levels(time)
    }
    class(fit) <- c("tpchreg", "pchreg", "phreg")
    fit
}

#' @export
extractAIC.tpchreg <- function(fit, scale, k = 2, ...){
    edf <- sum(fit$df) + length(fit$hazards)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

#' @export
nobs.tpchreg <- function(object, ...){
    object$n
}
