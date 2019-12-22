#' Proportional hazards regression with piecewise constant hazards and tabular data.
#' 
#' @param formula a formula with 'oe(count, exposure) ~ x1 + ...'
#' @param data a data frame with event, exposure, age plus covariates
#' @param time the time variable, a factor indicating time intervals.
#' @param pieces numeric vector of length 1 or length(levels(time)): The length(s)
#' of timeintervals.
#' @param subset subset of data, not implemented yet.
#' @param na.action Not implemented yet. 
#' @param contrasts Not implemented yet.
#' @param start.coef For themoment equal to zero.
#' @param control list of control parameters for the optimization. 

#' @note This function is under development and not well ducumented for the time
#' being. Use it with care, but it should work with standard (default) settings.
#' 
#' @seealso \code{\link{oe}}.
#' 
#' @export
tpchreg <- function(formula,
                    data,
                    time,
                    pieces, # length of pieces
                    subset,
                    na.action,
                    contrasts = NULL,
                    start.coef = NULL,
                    control = list(epsilon = 1.e-8,
                                   maxit = 200, trace = FALSE)){

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
    temp <- c("", "formula", "data", "time", "subset", "na.action")
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
    ##offset <- log(Y[, 2])
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

    time <- m$"(time)"
    if(is.null(strats)){
        strata <- rep(1, length(count))
    }else{
        strata <- strats
    }

    if (is.null(time)) time <- rep(1, length(count))

    n.ivls <- length(unique(time))
    if (missing(pieces)) {
        pieces <- rep(1, n.ivls)
    }else{
        if (length(pieces) == 1){
            pieces <- rep(pieces, n.ivls)
        }else{
            if (length(pieces) != n.ivls){
                stop("Interval count mismatch")
            }
        }
    }
    ##cbind(count, offset, X, strata, time)

    offset <- rep(0, length(count))
    ## Calling the work horse:
    fit <- tpchreg.fit(X, count, exposure, offset, strata, time, pieces)
    ##Terms
    fit$covars <- covars
    fit$terms <- Terms
    ##fit$newTerms <- newTerms
    stru <- struct(m, covars, dropx, X, exposure)
    fit$isI <- stru$isI
    fit$isF <- stru$isF
    fit$isO <- stru$isO
    fit$levels <- stru$levels
    fit$call <- call
    fit$w.means <- stru$w.means
    fit$ttr <- stru$ttr
    fit$pieces <- pieces
    class(fit) <- "tpchreg"
    fit
}

#'@export
extractAIC.tpchreg <- function(fit, scale, k = 2, ...){
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}
