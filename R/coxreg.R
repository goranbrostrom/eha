#' Cox regression
#' 
#' Performs Cox regression with some special attractions, especially
#' \emph{sampling of risksets} and \emph{the weird bootstrap}.
#' 
#' The default method, \code{efron}, and the alternative, \code{breslow}, are
#' both the same as in \code{\link[survival]{coxph}} in package
#' \code{survival}. The methods \code{mppl} and \code{ml} are maximum
#' likelihood, discrete-model, based.
#' 
#' @usage coxreg(formula = formula(data), data = parent.frame(), weights,
#' subset, t.offset, na.action = getOption("na.action"), init = NULL, method =
#' c("efron", "breslow", "mppl", "ml"), control = list(eps = 1e-08, maxiter =
#' 25, trace = FALSE), singular.ok = TRUE, model = FALSE, center = NULL, x =
#' FALSE, y = TRUE, hazards = NULL, boot = FALSE, efrac = 0, geometric = FALSE,
#' rs = NULL, frailty = NULL, max.survs = NULL, coxph = TRUE)
#'
#' @param formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function.
#' @param data a data.frame in which to interpret the variables named in the
#' formula.
#' @param weights Case weights; time-fixed or time-varying.
#' @param subset An optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param t.offset Case offsets; time-varying.
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.  Default is
#' \code{options()$na.action}.
#' @param init vector of initial values of the iteration.  Default initial
#' value is zero for all variables.
#' @param method Method of treating ties, "efron" (default), "breslow", "mppl"
#' (maximum partial partial likelihood), or "ml" (maximum likelihood).
#' @param control a list with components \code{eps} (convergence criterion),
#' \code{maxiter} (maximum number of iterations), and \code{silent} (logical,
#' controlling amount of output). You can change any component without mention
#' the other(s).
#' @param singular.ok Not used
#' @param model Not used
#' @param center deprecated. See Details.
#' @param x Return the design matrix in the model object?
#' @param y return the response in the model object?
#' @param hazards deprecated. Was: Calculate baseline hazards? Default is TRUE.
#' Calculating hazards is better done separately, after fitting. In most cases.
#' @param boot Number of boot replicates. Defaults to FALSE, no boot samples.
#' @param efrac Upper limit of fraction failures in 'mppl'.
#' @param geometric If TRUE, forces an 'ml' model with constant riskset
#' probability. Default is FALSE.
#' @param rs Risk set?
#' @param frailty Grouping variable for frailty analysis. Not in use (yet).
#' @param max.survs Sampling of risk sets? If given, it should be (the upper
#' limit of) the number of survivors in each risk set.
#' @param coxph Logical, defaults to \code{TRUE}. Determines if standard work
#' should be passed to \code{\link[survival]{coxph}} via entry points. 
#' @return A list of class \code{c("coxreg", "coxph")} with components
#' \item{coefficients}{Fitted parameter estimates.}
#' \item{var}{Covariance matrix of the estimates.}
#' \item{loglik}{Vector of length two; first component is the value at
#' the initial parameter values, the second component
#' is the maximized value.}
#' \item{score}{The score test statistic (at the initial value).}
#' \item{linear.predictors}{The estimated linear predictors.}
#' \item{residuals}{The martingale residuals.}
#' \item{hazards}{The estimated baseline hazards, calculated at the value zero of
#' the covariates (rather, columns of the design matrix). Is a list,
#' with one component per stratum. Each
#' component is a matrix with two columns, the first contains risk times, the
#' second the corresponding hazard atom.}
#' \item{means}{Means of the columns of
#' the design matrix corresponding to covariates, if \code{center = TRUE}.
#' Columns corresponding to factor levels gice a zero in the corresponding
#' position in \code{means}. If \code{center = FALSE}, \code{means} are all
#' zero.}
#' \item{w.means}{Weighted (against exposure time) means of covariates;
#' weighted relative frequencies of levels of factors.}
#' \item{n}{Number of spells in indata (possibly after removal of cases
#' with NA's).}
#' \item{n.events}{Number of events in data.}
#' \item{terms}{Used by extractor functions.}
#' \item{assign}{Used by extractor functions.}
#' \item{y}{The Surv vector.}
#' \item{isF}{Logical vector indicating the covariates that are factors.}
#' \item{covars}{The covariates.}
#' \item{ttr}{Total Time at Risk.}
#' \item{levels}{List of levels of factors.}
#' \item{formula}{The calling formula.}
#' \item{bootstrap}{The (matrix of) bootstrap replicates, if requested on
#' input. It is up to the user to do
#' whatever desirable with this sample.}
#' \item{boot.sd}{The estimated standard errors of the bootstrap replicates.}
#' \item{call}{The call.}
#' \item{method}{The method.}
#' \item{convergence}{Did the optimization converge?}
#' \item{fail}{Did the optimization fail? (Is \code{NULL} if not).}
#' @note This function starts by creating risksets, if no riskset is supplied
#' via \code{rs}, with the aid of \code{\link{risksets}}. Supplying output from
#' \code{risksets} via \code{rs} fails if there are any NA's in the data! Note
#' also that it depends on stratification, so \code{rs} contains information
#' about stratification. Giving another strata variable in the formula is an
#' error. The same is ok, for instance to supply stratum interactions.
#' @section Warning: The use of \code{rs} is dangerous, see note. It can
#' however speed up computing time considerably for huge data sets.
#' @author Göran Broström
#' @seealso \code{\link[survival]{coxph}}, \code{\link{risksets}}
#' @references Broström, G. and Lindkvist, M. (2008). Partial partial
#' likelihood. Communications in Statistics: Simulation and Computation 37:4,
#' 679-686.
#' @keywords regression survival
#' @examples
#' 
#'  dat <- data.frame(time=  c(4, 3,1,1,2,2,3),
#'                 status=c(1,1,1,0,1,1,0),
#'                 x=     c(0, 2,1,1,1,0,0),
#'                 sex=   c(0, 0,0,0,1,1,1))
#'  coxreg( Surv(time, status) ~ x + strata(sex), data = dat) #stratified model
#'  # Same as:
#'  rs <- risksets(Surv(dat$time, dat$status), strata = dat$sex)
#'  coxreg( Surv(time, status) ~ x, data = dat, rs = rs) #stratified model
#'  
#' @export
coxreg <- function (formula = formula(data),
                    data = parent.frame(),
                    weights,
                    subset,
                    t.offset,
                    na.action = getOption("na.action"),
                    init = NULL,
                    method = c("efron", "breslow", "mppl", "ml"),
                    control = list(eps = 1e-8, maxiter = 25, trace = FALSE),
                    singular.ok = TRUE,
                    model = FALSE,
                    center = NULL,
                    x = FALSE,
                    y = TRUE,
                    hazards = NULL,
                    boot = FALSE,
                    efrac = 0,
                    geometric = FALSE,
                    rs = NULL,
                    frailty = NULL,
                    max.survs = NULL,
                    coxph = TRUE)
{

    if (!missing(center)){
        warning("argument 'center' is deprecated.")
    }
    if (!missing(hazards)){
        warning("argument 'hazards' is deprecated.")
    }
    
    meth <- method[1]
    if (coxph){
        cox.ph <- (missing(t.offset) &&
                       (meth %in% c("breslow", "efron")) &&
                       is.null(rs) &&
                       is.null(max.survs) &&
                       (!boot) &&
                       (efrac == 0) &&
                       is.null(frailty) &&
                       (!geometric))
        if (!cox.ph) warning("'coxph is not called despite 'coxph = TRUE'")
    }else{
        cox.ph <- coxph # == TRUE
    }
    if (FALSE){  ############################## NOTE!!!!###########
    ##if (cox.ph){
        Call <- match.call()
        Call[[1]] <- quote(survival::coxph)
        ##return(Call)
        fit <- eval.parent(Call)
        class(fit) <- c("coxreg", "coxph")
        ##return(fit)
    }
    if (!is.null(frailty))
        stop("Frailty not implemented (yet). Try the 'coxme' package")
    method <- match.arg(method)
    
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "weights", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]

    special <- "strata"
    Terms <- if (missing(data))
        terms(formula, special)
    else
        terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")

    if (is.null(max.survs)) max.survs <- NROW(Y)
    if (missing(weights)) weights <- rep(1, NROW(Y))
    else weights <- model.extract(m, "weights")
    cox.ph <- cox.ph && (length(weights) == NROW(Y))
    
    if (missing(t.offset)) t.offset <- NULL
    ##
    offset <- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0)
        rep(0, nrow(Y))
    else if (tt == 1)
        m[[offset]]
    else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
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
                     function(x) x -1)
    X <- X[, -1, drop = FALSE]
    ##}
#########################################
        
    if (length(dropx)){
        covars <- names(m)[-c(1, (dropx + 1))]
    }else{
        covars <- names(m)[-1]
    }

    isI <- logical(NCOL(X)) # Added Jan 2014; 2.4-0.
    isF <- logical(length(covars))
    isO <- logical(length(covars))
    if (length(covars)){
        for (i in 1:length(covars)){
            if (length(dropx)){
                if (is.logical(m[, -(dropx + 1)][, (i + 1)])){
                    m[, -(dropx + 1)][, (i + 1)] <-
                        as.factor(m[, -(dropx + 1)][, (i + 1)])
                }
                isF[i] <- is.factor(m[, -(dropx + 1)][, (i + 1)])## ||
                isO[i] <- is.ordered(m[, -(dropx + 1)][, (i + 1)])## ||
                ##is.logical(m[, -(dropx + 1)][, (i + 1)]) )
            }else{
                if (is.logical(m[, (i + 1)])){
                    m[, (i + 1)] <- as.factor(m[, (i + 1)])
                }
                isF[i] <- is.factor(m[, (i + 1)]) ##||
                isO[i] <- is.ordered(m[, (i + 1)]) ##||
                ## is.logical(m[, (i + 1)]) )
            }
        }
    }
    
    if (any(isF)){
        levels <- list()
        index <- 0
        for ( i in 1:length(covars) ){
            if (isF[i]){
                index <- index + 1
                if (length(dropx)){
                    ll <- levels(m[, -(dropx + 1)][, (i + 1)])
                    if (isO[i]){
                        levels[[i]] <- paste("order", seq_along(ll) - 1)
                    }else{
                        levels[[i]] <- ll
                    }
                }else{
                    ll <- levels(m[, (i + 1)])
                    if (isO[i]){
                        levels[[i]] <- paste("order", seq_along(ll) - 1)
                    }else{
                        levels[[i]] <- ll
                    }
                }
            }else{
                levels[[i]] <- NULL
            }
        }
    }else{
        levels <- NULL
    }

   if (any(isF)){ ## New; get isI: (Jan 2014; 2.4-0)
        indx <- 0
        for (i in 1:length(covars)){
            indx <- indx + 1
            if (isF[i]){
                isI[indx] <- TRUE
                if (length(levels[[i]]) >= 3){
                    for (j in 3:length(levels[[i]])){
                        indx <- indx + 1
                        isI[indx] <- TRUE
                    }
                }
            }
        }
    }

    n.events <- sum(Y[, NCOL(Y)] != 0)
    if (n.events == 0) stop("No events; no sense in continuing!")
    
    if (FALSE){   ########################## NOTE!!!! ################### 
    ##if (cox.ph){  ## NEW ++++
        fit$df <- length(fit$coefficients)
        fit$isF <- isF
        fit$isI <- isI
        fit$covars <- covars
        fit$n.events <- n.events
        fit$nullModel <- (NCOL(X) == 0)
        fit$levels <- levels
        fit$n.events <- n.events
        fit <- wMeans(fit, Y, m, isF)
        class(fit) <- "coxreg"
        if (length(strats)){
            fit$stratum <- levels(as.factor(strata.keep)) ## New
            ## 'stratum' to be out of the way for 'strata' in survfit!
        }
        return(fit)
    }
    ##########################################

    ## Fixed now? if (FALSE){      ## This has to be fixed in the future!!
    if (NCOL(X) == 0){ # No covariates; special treatment!
        if (is.null(strats)){
            stratum <- rep(1, NROW(Y))
        }else{
            stratum <- strata.keep
        }
        type <- attr(Y, "type")
        control$iter.max <- 0
        control$toler.chol <- .Machine$double.eps^0.75
        control$toler.inf <- sqrt(control$eps)
        control$outer.max <- 10
        X <- matrix(0, nrow = NROW(Y), ncol = 1)
        init <- 0
        if (FALSE){
        if (type == "counting"){
            fit <- survival::agreg.fit(X, Y, stratum, offset, init,
                                        control, weights = weights,
                                        method = method, row.names(m))
        }else{
            fit <- survival::coxph.fit(X, Y, stratum, offset, init,
                                        control, weights = weights,
                                        method = method, row.names(m))
        }
        } # End if (FALSE)
        fit$nullModel <- TRUE

        ##if (hazards){
        scores <- exp(offset)
        hazards <- getHaz(Y, stratum, scores) ## Make this fast with C?
        class(hazards) <- "hazdata"
        fit$hazards <- hazards
        ##}
        fit$call <- call
        fit$n.events <- n.events
        fit$n <- NROW(Y)
        fit$y <- Y
        class(fit) <- c("coxreg", "coxph")
        fit$means <- 0
        fit$terms <- Terms
        return(fit)
    
    }else if (cox.ph){
        type <- attr(Y, "type")
        control$iter.max <- control$maxiter
        control$toler.chol <- .Machine$double.eps^0.75
        control$toler.inf <- sqrt(control$eps)
        control$outer.max <- 10
        if (is.null(strats)) stratum <- rep(1, NROW(Y))
        else stratum <- strats
        if (type == "counting"){
            fit <- survival::agreg.fit(X, Y, stratum, offset, init,
                                        control, weights = weights,
                                        method = method, row.names(m))
        }else{
            fit <- survival::coxph.fit(X, Y, stratum, offset, init,
                                        control, weights = weights,
                                        method = method, row.names(m))
        }
        fit$nullModel <- FALSE
        ## get hazards
        ## New in 2.4-0: covariates are centered; indicators not!
        ## If center == FALSE, X.means are "added back" before call to
        ## getHaz!!

        if (FALSE){ # 'center' is deprecated.
        ##if (center){
            X.means <- colMeans(X)
            for (i in seq_len(NCOL(X))){
                if (isI[i]) X.means[i] <- 0
            }
            scores <- exp(offset + X %*% fit$coefficients -
                          sum(X.means * fit$coefficients))
        }else{
            X.means <- numeric(NCOL(X))
            scores <- exp(offset + X %*% fit$coefficients)
        }

        ##if (hazards){
         ##   hazards <- getHaz(Y, stratum, scores)
        ##    class(hazards) <- "hazdata"
        ##    fit$hazards <- hazards
        ##}else{
        ##    fit$hazards <- NULL
        ##}
            ##rs <- risksets(Y, strats)
        ##hazard <- .Fortran("gethaz",
        ##                   as.integer(NROW(Y)),  # 'nn'
        ##                   as.integer(length(rs$antrs)), # 'ns' 
        ##                   as.integer(rs$antrs), # 'antrs'
        ##                   as.integer(rs$size), # 'size'
        ##                   as.integer(rs$n.events), # 'nevents'
        ##                   as.integer(length(rs$riskset)), # 'totsize'
        ##                   as.integer(rs$riskset), # 'riskset'
        ##                   as.double(exp(fit$linear.predictors)), # 'score'
        ##                   as.integer(sum(rs$antrs)), # 'totrs'
        ##                   hazard = double(sum(rs$antrs)), # 'hazard' (return)
        ##                   DUP = FALSE,
        ##                   PACKAGE = "eha")$hazard
        ## Put it on:
        ##haz.mean <- fit$hazard::: At means of covariates:
        ##if (!is.null(fit$coefficients))
          ##  hazard <- 1 - (1 - hazard)^exp(fit$means * fit$coefficients)
        ##hazards <- list()
        ##stopp <- cumsum(rs$antrs)
        ##startt <- c(1, 1 + stopp[-length(rs$antrs)])
        ##for (i in 1:length(rs$antrs)){
          ##  hazards[[i]] <- cbind(rs$risktimes[startt[i]:stopp[i]],
            ##                      hazard[startt[i]:stopp[i]])
        ##}
        ##fit$hazards <- hazard
    }else{ # if (!cox.ph)
        if (NCOL(Y) == 2){
            Y <- cbind(numeric(NROW(Y)), Y)
            attr(Y, "type") <- "counting"
        }
        
        ##return(Y)
        type <- attr(Y, "type")
        if (type != "right" && type != "counting")
            stop(paste("Cox model doesn't support \"", type, "\" survival data",
                       sep = ""))
        
        if ((!is.null(init)) && (length(init) != NCOL(X)))
            stop("Wrong length of 'init'")
        
        
        if (is.list(control)){
            if (is.null(control$eps)) control$eps <- 1e-8
            if (is.null(control$maxiter)) control$maxiter <- 10
            if (is.null(control$trace)) control$trace <- FALSE
        }else{
            stop("control must be a list")
        }
        
### New start for cox.ph (not any more!) ##################################
        if (geometric){
            method <- "ml"
            fit <- geome.fit(X,
                             Y,
                             rs,
                             strats,
                             offset,
                             init,
                             max.survs,
                             method,
                             ##                         boot,
                             control)
            fit$nullModel <- FALSE
        }else{
            fit <- coxreg.fit(X,
                              Y,
                              rs,
                              weights,
                              t.offset,
                              strats,
                              offset,
                              init,
                              max.survs,
                              method,
                              center = NULL, # Deprecated
                              boot,
                              efrac,
                              calc.hazards = NULL, # Deprecated
                              calc.martres = TRUE,
                              control,
                              verbose = TRUE)
        }
        ## get hazards
        fit$nullModel <- FALSE
        if (FALSE){ # 'center' is deprecated
        ##if (center){
            X.means <- colMeans(X)
            for (i in seq_len(NCOL(X))){
                if (isI[i]) X.means[i] <- 0
            }
            scores <- exp(offset + X %*% fit$coefficients -
                          sum(X.means * fit$coefficients))
        }else{
            X.means <- numeric(NCOL(X))
            scores <- exp(offset + X %*% fit$coefficients)
        }
        
        if (is.null(strats)){
            stratum <- rep(1, NROW(Y))
        }else{
            stratum <- strats
            fit$stratum <- strats
        }
        ##if (hazards){
        ##    hazards <- getHaz(Y, stratum, scores)
        ##    class(hazards) <- "hazdata"
        ##    fit$hazards <- hazards
        ##}else{
        ##    fit$hazards <- NULL
        ##}
        

        fit$convergence <- as.logical(fit$conver)
        fit$conver <- NULL ## Ugly!
        fit$f.convergence <- as.logical(fit$f.conver)
        fit$f.conver <- NULL
    }
###########################################################################
## Crap dealt with ......

    if (is.character(fit)) {
        fit <- list(fail = fit)
        class(fit) <- "coxreg"
    }else if ((!cox.ph) && !fit$fail){
        if (length(fit$coef) && any(is.na(fit$coef))) {
            vars <- (1:length(fit$coef))[is.na(fit$coef)]
            msg <- paste("X matrix deemed to be singular; variable",
                         paste(vars, collapse = " "))
            if (singular.ok)
                warning(msg)
            else stop(msg)
        }
        fit$n <- nrow(Y)
        class(fit) <- fit$method
        fit$terms <- Terms
        fit$assign <- assign
        if (FALSE){ ## Out-commented
            if (length(fit$coef) && is.null(fit$wald.test)) {
                nabeta <- !is.na(fit$coef)
                if (is.null(init))
                    temp <- fit$coef[nabeta]
                else temp <- (fit$coef - init)[nabeta]
                ##fit$wald.test <-
                  ##  survival:::coxph.wtest(fit$var[nabeta, nabeta],
                    ##                       temp, control$toler.chol)$test
            }
        } ## End Out-commented
        na.action <- attr(m, "na.action")
        if (length(na.action))
            fit$na.action <- na.action
        if (model)
            fit$model <- m
        if (x) {
            fit$x <- X
            if (length(strats))
                fit$strata <- strata.keep
        }
    }
    fit$stratum <- strats
    if (y)
        fit$y <- Y
    if (x)
        fit$x <- X
    ##if (!is.null(weights) && any(weights != 1))
    ##    fit$weights <- weights

    ##########################################

    fit$isI <- isI
    fit$isF <- isF
    fit$isO <- isO
    fit$covars <- covars
    fit$levels <- levels
    fit <- wMeans(fit, Y, m, isF)
    
    ##########################################
    fit$nullModel <- FALSE

    fit$formula <- formula(Terms)
    fit$terms <- Terms
    fit$call <- call
    fit$n.events <- n.events
    ##fit$center <- center
    if (length(fit$coefficients)){
        names(fit$coefficients) <- colnames(X)
        fit$means <- X.means
    }else{
        fit$means <- numeric(0)
    }
    if (length(strats)){
        fit$strata <- levels(as.factor(strata.keep)) ## New 2.2-6
    }
    fit$method <- method
    fit$n <- NROW(Y)
    fit$df <- length(fit$coefficients)
    fit$lin.pred <- fit$linear.predictors # Just for error check ...
    fit$linear.predictors <- offset + X %*% fit$coefficients
    
    class(fit) <- c("coxreg", "coxph") # Not Removed "coxph"; cox.zph!
    ##class(fit) <- "coxreg"
    fit
}

wMeans <- function(fit, Y, m, isF){
    if (NCOL(Y) == 3){
        s.wght <- (Y[, 2] - Y[, 1])## * weights
    }else{
        s.wght <- Y[, 1]
    }
    fit$ttr <- sum(s.wght)
    fit$w.means <- list()
    if (length(fit$covars)){
        for (i in 1:length(fit$covars)){
            nam <- fit$covars[i]
            col.m <- which(nam == names(m))
            if (isF[i]){
                n.lev <- length(fit$levels[[i]])
                fit$w.means[[i]] <- numeric(n.lev)
                for (j in 1:n.lev){
                    who <- m[, col.m] == fit$levels[[i]][j]
                    fit$w.means[[i]][j] <-
                        sum( s.wght[who] ) / fit$ttr ## * 100, if in per cent
                }
            }else{
                fit$w.means[[i]] <- sum(s.wght * m[, col.m]) / fit$ttr
            }
        }
    }
    fit
}
