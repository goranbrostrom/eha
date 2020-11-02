#' Weibull Regression
#' 
#' Proportional hazards model with baseline hazard(s) from the Weibull family
#' of distributions.  Allows for stratification with different scale and shape
#' in each stratum, and left truncated and right censored data.
#' 
#' The parameterization is the same as in \code{\link{coxreg}} and
#' \code{\link[survival]{coxph}}, but different from the one used by
#' \code{\link[survival]{survreg}}. The model is \deqn{h(t; a, b, \beta, z) =
#' (a/b) (t/b)^{a-1} exp(z\beta)}{% h(t; a, b, beta, z) = (a/b) (t/b)^(a-1)
#' exp(z beta)} This is in correspondence with \code{\link{Weibull}}. To
#' compare regression coefficients with those from \code{survreg} you need to
#' divide by estimated shape (\eqn{\hat{a}}{a}) and change sign. The p-values
#' and test statistics are however the same, with one exception; the score test
#' is done at maximized scale and shape in \code{weibreg}.
#' 
#' This model is a Weibull distribution with shape parameter \eqn{a} and scale
#' parameter \eqn{b \exp(-z\beta / a)}{b exp(-z beta / a)}
#' 
#' @param formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right.  The response must be a survival
#' object as returned by the Surv function.
#' @param data a data.frame in which to interpret the variables named in the
#' formula.
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.  Default is
#' \code{options()$na.action}.
#' @param init vector of initial values of the iteration.  Default initial
#' value is zero for all variables.
#' @param shape If positive, the shape parameter is fixed at that value (in
#' each stratum).  If zero or negative, the shape parameter is estimated.  If
#' more than one stratum is present in data, each stratum gets its own
#' estimate.
#' @param control a list with components \code{eps} (convergence criterion),
#' \code{maxiter} (maximum number of iterations), and \code{silent} (logical,
#' controlling amount of output).  You can change any component without mention
#' the other(s).
#' @param singular.ok Not used.
#' @param model Not used.
#' @param x Return the design matrix in the model object?
#' @param y Return the response in the model object?
#' @param center Deprecated, and not used.  Will be removed in the future.
#' @return A list of class \code{c("weibreg", "coxreg")} with components
#' \item{coefficients}{Fitted parameter estimates.} \item{var}{Covariance
#' matrix of the estimates.} \item{loglik}{Vector of length two; first
#' component is the value at the initial parameter values, the second componet
#' is the maximized value.} \item{score}{The score test statistic (at the
#' initial value).} \item{linear.predictors}{The estimated linear predictors.}
#' \item{means}{Means of the columns of the design matrix.}
#' \item{w.means}{Weighted (against exposure time) means of covariates;
#' weighted relative frequencies of levels of factors.} \item{n}{Number of
#' spells in indata (possibly after removal of cases with NA's).}
#' \item{events}{Number of events in data.} \item{terms}{Used by extractor
#' functions.} \item{assign}{Used by extractor functions.} %
#' \item{wald.test}{The Wald test statistic (at the initial value).}
#' \item{y}{The Surv vector.} \item{isF}{Logical vector indicating the
#' covariates that are factors.} \item{covars}{The covariates.}
#' \item{ttr}{Total Time at Risk.} \item{levels}{List of levels of factors.}
#' \item{formula}{The calling formula.} \item{call}{The call.}
#' \item{method}{The method.} \item{convergence}{Did the optimization
#' converge?} \item{fail}{Did the optimization fail? (Is \code{NULL} if not).}
#' \item{pfixed}{TRUE if shape was fixed in the estimation.}
#' @note This function is not maintained, and may behave in unpredictable ways.
#' Use \code{\link{phreg}} with \code{dist = "weibull"} (the default) instead!
#' Will soon be declared deprecated.
#' @section Warning: The print method \code{\link{print.weibreg}} doesn't work
#' if threeway or higher order interactions are present.
#' 
#' Note further that covariates are internally centered, if \code{center =
#' TRUE}, by this function, and this is not corrected for in the output. This
#' affects the estimate of \eqn{\log(scale)}{log(scale)}, but nothing else. If
#' you don't like this, set \code{center = FALSE}.
#' @author Göran Broström
#' @seealso \code{\link{phreg}}, \code{\link{coxreg}},
#' \code{\link{print.weibreg}}
#' @keywords survival regression
#' @examples
#' 
#'  dat <- data.frame(time = c(4, 3, 1, 1, 2, 2, 3),
#'                 status = c(1, 1, 1, 0, 1, 1, 0),
#'                 x = c(0, 2, 1, 1, 1, 0, 0),
#'                 sex = c(0, 0, 0, 0, 1, 1, 1))
#'  weibreg( Surv(time, status) ~ x + strata(sex), data = dat) #stratified model
#' 
#' @export weibreg
weibreg <-
  function (formula = formula(data),
            data = parent.frame(), 
            na.action = getOption("na.action"),
            init,
            shape = 0, ## Means shape is estimated, ie true Weibull; > 0 fixed!
            control = list(eps = 1e-4, maxiter = 10, trace = FALSE),
            singular.ok = TRUE,
            model = FALSE, 
            x = FALSE,
            y = TRUE,
            center = TRUE) 
{
    
    pfixed <- (shape > 0)
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    
    special <- "strata"
    Terms <- if (missing(data)) 
      terms(formula, special)
    else terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    
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
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 1)
    X <- X[, -1, drop = FALSE]
    ncov <- NCOL(X)
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
                isF[i] <- ( is.factor(m[, -(dropx + 1)][, (i + 1)]) ||
                           is.logical(m[, -(dropx + 1)][, (i + 1)]) )
            }else{
                isF[i] <- ( is.factor(m[, (i + 1)]) ||
                           is.logical(m[, (i + 1)]) )
            }      
        }
        
        if (ant.fak <- sum(isF)){
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
                    ##cat("NULL level no  ", i, "\n")
                    levels[[i]] <- NULL
                }
            }
        }else{
            levels <- NULL
        }
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
    if (missing(init)) 
      init <- NULL
    
    if (is.list(control)){
        if (is.null(control$eps)) control$eps <- 1e-4
        if (is.null(control$maxiter)) control$maxiter <- 10
        if (is.null(control$trace)) control$trace <- FALSE
    }else{
        stop("control must be a list")
    }
    
    
    fit <- weibreg.fit(X, 
                       Y,
                       strats,
                       offset,
                       init,
                       shape,
                       control,
                       center)

    if (ncov){
        fit$linear.predictors <- offset + X %*% fit$coefficients[1:ncov]
        fit$means <- apply(X, 2, mean)
    }else{
        fit$linear.predictors <- NULL
        fit$means <- NULL
    }
    ##score <- exp(lp)


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
                fit$w.means[[i]] <- sum(s.wght * m[, col.m]) / fit$ttr
            }
        }
        fit$means <- colMeans(X)
        
    }
    
##########################################
    fit$ttr <- sum(s.wght)
    ##names(fit$coefficients) <- coef.names 
    fit$levels <- levels
    fit$formula <- formula(Terms)
    fit$call <- call
    fit$n.events <- n.events 
    ##class(fit) <- c("weibreg", "coxreg", "coxph")
    class(fit) <- c("weibreg", "phreg")
    fit$pfixed <- pfixed
    fit
}
