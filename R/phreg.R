#' Parametric Proportional Hazards Regression
#' 
#' Proportional hazards model with parametric baseline hazard(s).  Allows for
#' stratification with different scale and shape in each stratum, and left
#' truncated and right censored data.
#' 
#' The parameterization is the same as in \code{\link{coxreg}} and
#' \code{\link[survival]{coxph}}, but different from the one used by
#' \code{\link[survival]{survreg}} (which is not a proportional hazards
#' modelling function). The model is \deqn{S(t; a, b, \beta, z) =
#' S_0((t/b)^a)^{\exp((z-mean(z))\beta)}}{% S(t; a, b, beta, z) =
#' S0((t/b)^a)^exp((z - mean(z)) beta)} where S0 is some standardized survivor
#' function.
#' 
#' @param formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right.  The response must be a survival
#' object as returned by the Surv function.
#' @param data a data.frame in which to interpret the variables named in the
#' formula.
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.  Default is
#' \code{options()$na.action}.
#' @param dist Which distribution? Default is "weibull", with the alternatives
#' "ev" (Extreme value), "gompertz", "pch" (piecewise constant hazards
#' function), "loglogistic" and "lognormal". A special case like the
#' \code{exponential} can be obtained by choosing "weibull" in combination with
#' \code{shape = 1}, or "pch" without \code{cuts}.
#' @param cuts Only used with \code{dist = "pch"}. Specifies the points in time
#' where the hazard function jumps. If omitted, an exponential model is fitted.
#' @param init vector of initial values of the iteration.  Default initial
#' value is zero for all variables.
#' @param shape If positive, the shape parameter is fixed at that value (in
#' each stratum).  If zero or negative, the shape parameter is estimated.  If
#' more than one stratum is present in data, each stratum gets its own
#' estimate. Only relevant for the Weibull and Extreme Value distributions.
#' @param param Applies only to the Gompertz distribution: "canonical" is
#' defined in the description of the \code{\link{Gompertz}} distribution;
#' "rate" transforms \code{scale} to 1/log(scale), giving the same
#' parametrization as in Stata and SAS. The latter thus allows for a negative
#' rate, or a "cure" (Gompertz) model. The default is "canonical"; if this
#' results in extremely large scale and/or shape estimates, consider trying
#' "rate".
#' @param control a list with components \code{eps} (convergence criterion),
#' \code{maxiter} (maximum number of iterations), and \code{silent} (logical,
#' controlling amount of output).  You can change any component without mention
#' the other(s).
#' @param singular.ok Not used.
#' @param model Not used.
#' @param x Return the design matrix in the model object?
#' @param y Return the response in the model object?
#' @param center Is now deprecated. 
#' @return A list of class \code{c("phreg", "coxreg")} with components
#' \item{coefficients}{Fitted parameter estimates.} 
#' \item{cuts}{Cut points for
#' the "pch" distribution. \code{NULL} otherwise.} 
#' \item{hazards}{The estimated
#' constant levels in the case of the "pch" distribution. \code{NULL}
#' otherwise.} 
#' \item{var}{Covariance matrix of the estimates.}
#' \item{loglik}{Vector of length two; first component is the value at the
#' initial parameter values, the second componet is the maximized value.}
#' \item{score}{The score test statistic (at the initial value).}
#' \item{linear.predictors}{The estimated linear predictors.}
#' \item{means}{Means of the columns of the design matrix, except those columns
#' corresponding to a factor level. Otherwise all
#' zero.} 
#' \item{w.means}{Weighted (against exposure time) means of covariates;
#' weighted relative frequencies of levels of factors.} 
#' \item{n}{Number of
#' spells in indata (possibly after removal of cases with NA's).}
#' \item{n.events}{Number of events in data.} 
#' \item{terms}{Used by extractor functions.} 
#' \item{assign}{Used by extractor functions.} %
#' \item{wald.test}{The Wald test statistic (at the initial value).}
#' \item{y}{The Surv vector.} 
#' \item{isF}{Logical vector indicating the
#' covariates that are factors.} 
#' \item{covars}{The covariates.}
#' \item{ttr}{Total Time at Risk.} 
#' \item{levels}{List of levels of factors.}
#' \item{formula}{The calling formula.} 
#' \item{call}{The call.}
#' \item{method}{The method.} 
#' \item{convergence}{Did the optimization
#' converge?} 
#' \item{fail}{Did the optimization fail? (Is \code{NULL} if not).}
#' \item{pfixed}{TRUE if shape was fixed in the estimation.}
#' @note The lognormal and loglogistic baseline distributions are extended to a
#' three-parameter family by adding a "proportionality" parameter (multiplying
#' the baseline hazard function). The log of the estimated parameter turns up
#' as '(Intercept)' in the printed output. The reason for this extension is
#' that the standard lognormal and loglogistic distributions are not closed
#' under proportional hazards.
#' @section Warning: The lognormal and loglogistic distributions are included
#' on an experimental basis for the moment. Use with care, results may be
#' unreliable!
#' 
#' The gompertz distribution has an exponentially increasing hazard function
#' under the canonical parametrization. This may cause instability in the
#' convergence of the fitting algorithm in the case of near-exponential data.
#' It may be resolved by using \code{param = "rate"}.
#' @author Göran Broström
#' @seealso \code{\link{coxreg}}, \code{\link{check.dist}},
#' \code{link{aftreg}}.
#' @keywords survival regression
#' @examples
#' 
#' data(mort)
#' fit <- phreg(Surv(enter, exit, event) ~ ses, data = mort)
#' fit
#' plot(fit)
#' fit.cr <- coxreg(Surv(enter, exit, event) ~ ses, data = mort)
#' check.dist(fit.cr, fit)
#' 
#' @export phreg
phreg <- function (formula = formula(data),
                   data = parent.frame(),
                   na.action = getOption("na.action"),
                   dist = "weibull",
                   cuts = NULL,
                   init,
                   shape = 0,
                   ## Means shape is estimated, ie true Weibull; > 0 fixed!
                   param = c("canonical", "rate"),
                   control = list(eps = 1e-8, maxiter = 20, trace = FALSE),
                   singular.ok = TRUE,
                   model = FALSE,
                   x = FALSE,
                   y = TRUE,
                   center = NULL) # NOTE: Changed from 'NULL' in 1.4-1
{                                 # NOTE: Changed back in 2.2-2!
                                 ## NOTE: Changed again in 2.4-0; affects only 
                                 ## plot and log(scale)
                                 ## Finally (2.8.2) deprecated. 
                                 ## Centering never a good idea!
    if (!missing(center)){
      warning("argument 'center' is deprecated: 
              Reported results are not centered.", call. = FALSE)
    }
    param <- param[1]
    pfixed <- any(shape > 0) & dist %in% c("weibull", "ev") # Fixed 2020-07-26!
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)] # m is a call

    special <- "strata"
    Terms <- if (missing(data))
      terms(formula, special)
    else terms(formula, special, data = data) # Terms is a 'terms' 'formula'

    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    ##return(m)
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv"))
      stop("Response must be a survival object")

    weights <- model.extract(m, "weights")
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
        ##if (dist == "pch") # Changed 2.4-0
          ##  stop("No strata allowed in the pch model (yet)") 
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
    ##return(X)
    if (!(dist %in% c("lognormal", "loglogistic"))){
        assign <- lapply(survival::attrassign(X, newTerms)[-1], function(x) x - 1)
    }else{
        assign <- lapply(survival::attrassign(X, newTerms), function(x) x - 1)
    }
    ## Lognormal & loglogistic need an intercept:

    if (!(dist %in% c("lognormal", "loglogistic"))){
        X <- X[, -1, drop = FALSE]
        intercept <- FALSE
        ncov <- NCOL(X)
    }else{
        intercept <- TRUE
        ncov <- NCOL(X) - 1
    }
    nullModel <- ncov == 0
#########################################

    if (ncov){
        if (length(dropx)){
            covars <- names(m)[-c(1, (dropx + 1))]
         }else{
             covars <- names(m)[-1]
        }

        isF <- logical(length(covars))
        for (i in 1:length(covars)){
            if (length(dropx)){
                if (is.logical(m[, -(dropx + 1)][, (i + 1)])){
                    m[, -(dropx + 1)][, (i + 1)] <-
                        as.factor(m[, -(dropx + 1)][, (i + 1)])
                }

                isF[i] <- is.factor(m[, -(dropx + 1)][, (i + 1)])## ||
                           ##is.logical(m[, -(dropx + 1)][, (i + 1)]) )
            }else{
                if (is.logical(m[, (i + 1)])){
                    m[, (i + 1)] <- as.factor(m[, (i + 1)])
                }
                
                isF[i] <- is.factor(m[, (i + 1)])## ||
                           ##is.logical(m[, (i + 1)]) )
            }
        }
         
        ant.fak <- sum(isF)
    }else{ #!ncov
        isF <- logical(0)
        ant.fak <- 0
    }
    ##cat("ant.fak = ", ant.fak, "\n")
    if (ant.fak){
        levels <- list()
        index <- 0
        for ( i in 1:length(covars) ){
            if (isF[i]){
                index <- index + 1
                if (length(dropx)){
                    levels[[i]] <- levels(m[, -(dropx + 1)][, (i + 1)])
                }else{
                    levels[[i]] <- levels(m[, (i + 1)])
                }
            }else{
                levels[[i]] <- NULL
            }
        }
    }else{
        levels <- NULL
    }

    isI <- logical(NCOL(X))
    if (ant.fak){
        indx <- 0
        for (i in seq_len(length(covars))){
            indx <- indx + 1
            if (isF[i]){
                isI[indx] <- TRUE
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
    if (FALSE){
    ##if (center){
        X.means <- colMeans(X)
        X.means[isI] <- 0
    }else{
        X.means <- 0
    }
    
##########################################
    type <- attr(Y, "type")
    if (type != "right" && type != "counting")
      stop(paste("This model doesn't support \"", type, "\" survival data",
                 sep = ""))

    if (NCOL(Y) == 2){
        Y <- cbind(numeric(NROW(Y)), Y)
    }

    n.events <- sum(Y[, 3] != 0)
    if (n.events == 0) stop("No events; no sense in continuing!")
    if (missing(init)){ # Is this wise? No!
        if (ncov){
            ##init <- coxreg(formula, data = data)$coefficients
            init <- rep(0, ncov)
        }else{
             init <- numeric(0)
         }
    }
    if (is.list(control)){
        if (is.null(control$eps)) control$eps <- 1e-8
        if (is.null(control$maxiter)) control$maxiter <- 10
        if (is.null(control$trace)) control$trace <- FALSE
    }else{
        stop("control must be a list")
    }
    if (dist == "gompertz"){
        if (param == "canonical"){
            
            fit <- gompreg(X,
                           Y,
                           strats,
                           offset,
                           init,
                           control)
        }else if (param == "rate"){
            fit <- gompregRate(X,
                               Y,
                               strats,
                               offset,
                               init,
                               control)
        }else{
            stop(paste(param, " is not a known parametrization."))
        }
        
        }else if(dist == "pch"){
            if (missing(cuts)){
                ##stop("'dist = pch' needs 'cuts' to be set")
                cuts <- numeric(0) # Exponential distribution(s)
            }
##        fit <- pchreg(X,
  ##                    Y,
    ##                  cuts,
      ##                offset,
        ##              init,
          ##            control,
            ##          center)
            fit <- pchreg2(X,
                           Y,
                           cuts,
                           offset,
                           strats,
                           init,
                           control,
                           center)
    }else{
        fit <- phreg.fit(X,
                         Y,
                         dist,
                         strats,
                         offset,
                         init,
                         shape,
                         control,
                         center = NULL)
    }
    
    if (fit$fail){
        warning(paste("Failed with error code ", fit$fail))
        return(1)
    }

    if (ncov){
      ##cat("ncov = ", ncov, "\n")
      fit$linear.predictors <- offset + X %*%
        fit$coefficients[1:(ncov + intercept)]
      ##fit$means <- X.means
    }else{
      fit$linear.predictors <- numeric(0)
      ##fit$means <- numeric(0)
    }
    ##score <- exp(lp)

    ##cat("fit$means == ", fit$means, "\n")
    if (!fit$fail){
        fit$fail <- NULL
    }else{
        out <- paste("Singular hessian; suspicious variable No. ",
                     as.character(fit$fail), ":\n",
                     names(coefficients)[fit$fail], " = ",
                     as.character(fit$value),
                     "\nTry running with fixed shape", sep = "")
        stop(out)
    }


    fit$convergence <- as.logical(fit$conver)
    fit$conver <- NULL ## Ugly!

###########################################################################
    ## Crap dealt with ......

    if (is.character(fit)) {
        fit <- list(fail = fit)
        class(fit) <- "mlreg"
    }
    else if (is.null(fit$fail)){
        if (!is.null(fit$coef) && any(is.na(fit$coef))) {
            vars <- (1:length(fit$coef))[is.na(fit$coef)]
            msg <- paste("X matrix deemed to be singular; variable",
                         paste(vars, collapse = " "))
            if (singular.ok)
              warning(msg)
            else stop(msg)
        }
        fit$n <- nrow(Y)
        fit$terms <- Terms
        fit$assign <- assign
        if (FALSE){ ## Out-commented...(why?)
            if (length(fit$coef) && is.null(fit$wald.test)) {
                nabeta <- !is.na(fit$coef)
                if (is.null(init))
                  temp <- fit$coef[nabeta]
                else temp <- (fit$coef - init)[nabeta]
                ##fit$wald.test <-
                  ##survival:::coxph.wtest(fit$var[nabeta, nabeta],
                    ##                     temp, control$toler.chol)$test
            }
        }
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
        if (y)
          fit$y <- Y
    }
    ##if (!is.null(weights) && any(weights != 1))
    ##    fit$weights <- weights

##########################################
    s.wght <- (Y[, 2] - Y[, 1])## * weights
    fit$ttr <- sum(s.wght)
    if (ncov){
        fit$isF <- isF
        fit$covars <- covars
        fit$w.means <- list()
        for (i in 1:length(fit$covars)){
            nam <- fit$covars[i]
            col.m <- which(nam == names(m))
            if (isF[i]){
                n.lev <- length(levels[[i]])
                fit$w.means[[i]] <- numeric(n.lev)
                for (j in 1:n.lev){
                    who <- m[, col.m] == levels[[i]][j]
                    fit$w.means[[i]][j] <-
                      sum( s.wght[who] ) / fit$ttr ## * 100, if in per cent
                }
            }else{
                ##fit$w.means[[i]] <- sum(s.wght * m[, col.m]) / fit$ttr
                fit$w.means[[i]] <- weighted.mean(m[, col.m], s.wght)
            }
        }
    }

##########################################
    fit$ttr <- sum(s.wght)
    ##names(fit$coefficients) <- coef.names
    fit$levels <- levels
    fit$formula <- formula(Terms)
    fit$call <- call
    fit$dist <- dist
    fit$n.events <- n.events
    fit$nullModel <- nullModel # Added 2020-07-26
    ##class(fit) <- c("phreg", "weibreg", "coxreg", "coxph")
    if (fit$dist == "pch"){
        class(fit) <- c("pchreg", "phreg")
    }else{
        class(fit) <- "phreg"
    }
    if (dist %in% c("weibull", "cv")){
        fit$pfixed <- pfixed
    }else{
        fit$pfixed <- FALSE
    }
    if (length(strats))
        fit$strata <- names(strats)
    if (length(strats)){
        fit$strata <- levels(as.factor(strata.keep))
    }
    fit
}
