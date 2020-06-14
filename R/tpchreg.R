#' Proportional hazards regression with piecewise constant hazards and tabular 
#' data.
#' 
#' @usage tpchreg(formula, data, time, subset, na.action, contrasts = NULL,
#' start.coef = NULL, 
#' control = list(epsilon = 1.e-8, maxit = 200, trace = FALSE))
#' 
#' @param formula a formula with 'oe(count, exposure) ~ x1 + ...'
#' @param data a data frame with event, exposure, age plus covariates
#' @param time the time variable, a factor indicating time intervals.
#' @param weights Case weights.
#' @param subset subset of data, not implemented yet.
#' @param na.action Not implemented yet.
#' @param contrasts Not implemented yet.
#' @param start.coef For the moment equal to zero.
#' @param control list of control parameters for the optimization. 

#' @note 
#' 
#' @seealso \code{\link{oe}}.
#' 
#' @export
tpchreg <- function(formula,
                    data,
                    time,
                    weights,
                    subset,
                    na.action,
                    contrasts = NULL,
                    start.coef = NULL,
                    control = list(epsilon = 1.e-8,
                                   maxit = 200, trace = FALSE)){
    if (missing(time)){
        stop("Argument 'time' is missing with no default")
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
        temp <- untangle.specials(Terms, "strata", 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars) == 1) 
            strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
        strats <- as.numeric(strata.keep)
    }
    if (length(dropx)) 
        newTerms <- Terms[-dropx]
    else newTerms <- Terms
    X <- model.matrix(newTerms, m)
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 
        1)
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
    
    time <- m$"(time)"
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

    ##if (is.null(time)) time <- rep(1, length(count))
    if (is.character(time)) {
        time <- as.factor(time)
    }
    if (!is.factor(time)){
        stop("Argument 'time' must be a character or factor variable")
    }
    cuts <- as.numeric(unique(unlist(strsplit(levels(time), "-"))))
    
    n.ivls <- length(unique(time))
   ## Calling the work horse:
    fit <- tpchreg.fit(X, count, exposure, offset, weights, stratum, time)
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
    fit$events <- sum(count)
    fit$ttr <- stru$ttr
    fit$n <- length(count)
    fit$cuts <- cuts
    fit$strata <- strata
    class(fit) <- c("tpchreg", "phreg")
    fit
}

#' @export
extractAIC.tpchreg <- function(fit, scale, k = 2, ...){
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

#' @export
nobs.tpchreg <- function(object, ...){
    object$n
}
