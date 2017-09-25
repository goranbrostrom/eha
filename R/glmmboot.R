#' Generalized Linear Models with fixed effects grouping
#' 
#' Fits grouped GLMs with fixed group effects. The significance of the grouping
#' is tested by simulation, with a bootstrap approach.
#' 
#' The simulation is performed by simulating new response vectors from the
#' fitted probabilities without clustering, and comparing the maximized log
#' likelihoods. The maximizations are performed by profiling out the grouping
#' factor. It is a very fast procedure, compared to \code{\link{glm}}, when the
#' grouping factor has many levels.
#' 
#' @param formula a symbolic description of the model to be fit. The details of
#' model specification are given below.
#' @param family Currently, the only valid values are \code{binomial} and
#' \code{poisson}. The binomial family allows for the \code{logit} and
#' \code{cloglog} links.
#' @param data an optional data frame containing the variables in the model.
#' By default the variables are taken from `environment(formula)', typically
#' the environment from which `glmmML' is called.
#' @param cluster Factor indicating which items are correlated.
#' @param weights Case weights.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action See glm.
#' @param offset this can be used to specify an a priori known component to be
#' included in the linear predictor during fitting.
#' @param start.coef starting values for the parameters in the linear
#' predictor.  Defaults to zero.
#' @param control Controls the convergence criteria. See
#' \code{\link{glm.control}} for details.
#' @param boot number of bootstrap replicates. If equal to zero, no test of
#' significance of the grouping factor is performed.
#' @return The return value is a list, an object of class 'glmmboot'.
#' \item{coefficients}{Estimated regression coefficients} \item{logLik}{the max
#' log likelihood} \item{cluster.null.deviance}{Deviance without the
#' clustering} \item{frail}{The estimated cluster effects} \item{bootLog}{The
#' logLik values from the bootstrap samples} \item{bootP}{Bootstrap p value}
#' \item{variance}{Variance covariance matrix} \item{sd}{Standard error of
#' regression parameters} \item{boot_rep}{No. of bootstrap replicates}
#' \item{mixed}{Logical} \item{deviance}{Deviance} \item{df.residual}{Its
#' degrees of freedom} \item{aic}{AIC} \item{boot}{Logical} \item{call}{The
#' function call}
#' @note There is no overall intercept for this model; each cluster has its own
#' intercept. See \code{frail}
#' @author Göran Broström
#' @seealso \code{link{glmmML}}, \code{\link{optim}}, \code{lmer} in the
#' package \code{lme4}, and \code{glmmPQL} in the package \code{MASS}.
#' @keywords regression nonlinear
#' @examples
#' 
#' ## Not run:
#' id <- factor(rep(1:20, rep(5, 20)))
#' y <- rbinom(100, prob = rep(runif(20), rep(5, 20)), size = 1)
#' x <- rnorm(100)
#' dat <- data.frame(y = y, x = x, id = id)
#' res <- glmmboot(y ~ x, cluster = id, data = dat, boot = 500)
#' ## End(Not run)
#' ##system.time(res.glm <- glm(y ~ x + id, family = binomial))
#' 
#' @export glmmboot
glmmboot <- function(formula,
                     family = binomial,
                     data,
                     cluster,
                     weights,
                     subset,
                     na.action,
                     offset,
                     start.coef = NULL,
                     control = list(epsilon = 1.e-8,
                     maxit = 200, trace = FALSE),
                     boot = 0){

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

    cl <- match.call()

    if (is.character(family)) 
        family <- get(family)
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
    
    if (missing(data))
        data <- environment(formula)
    
    mf <- match.call(expand.dots = FALSE)
    ## get a copy of the call; result: a list.
    
    mf$family <- mf$start.coef <- mf$start.sigma <- NULL
    mf$control <- mf$maxit <- mf$boot <- mf$conditional <- NULL
    mf$n.points <- mf$start.coef <- NULL
    mf[[1]] <- as.name("model.frame") # turn into a call to model.frame
    mf <- eval(mf, environment(formula)) # run model.frame
    
    ## Pick out the parts.
    mt <-  attr(mf, "terms")
    
    
    xvars <- as.character(attr(mt, "variables"))[-1]
    if ((yvar <- attr(mt, "response")) > 0) 
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(mf[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    
    p <- NCOL(X)
    
    Y <- model.response(mf, "numeric")
    offset <- model.offset(mf)
 
    cluster <- mf$"(cluster)"
    
    no.cluster <- (missing(cluster) || is.null(cluster) ||
                   (length(unique(cluster)) <= 1))
    if (no.cluster){
        warning("No (or constant) 'cluster'; consider using 'glm'")
        return(NULL)
    }    

    ## if (NCOL(Y) >  1) stop("Response must be univariate")
    
    if (!is.null(offset) && length(offset) != NROW(Y)) 
        stop(paste("Number of offsets is", length(offset), ", should equal", 
                   NROW(Y), "(number of observations)"))

    if (missing(weights)) weights <- rep.int(1, NROW(Y))
    if (any(weights < 0)) stop("negative weights not allowed")
    
    ## Remove eventual intercept from X.
    ## Taken care of thru separate intercepts for each 'cluster'.
    ## Weck: NOW,  31 jan 2005, we change that! ##
    ## Weck: NOW,  6 jun 2006, we move it to glmmbootFit! ##
    
    ##if (!is.na(coli <- match("(Intercept)", colnames(X))))
    ##    X <- X[, -coli, drop = FALSE]

    res <- glmmbootFit(X, Y,
                       weights,
                       start.coef,
                       cluster,
                       offset,
                       family,
                       control,
                       boot)
    
    res$mixed <- FALSE # ??????????????????

    if (family$family == "binomial"){
        res$deviance <- 2 * (sum(Y * log(ifelse(Y == 0, 1, Y/res$fitted)) +
                                (1-Y) *  log(ifelse(Y == 1, 1,
                                             (1-Y) / (1 - res$fitted)))))
    }else{
        res$deviance <- 2 * sum(Y * log(ifelse(Y == 0, 1, Y / res$fitted)) -
                             (Y - res$fitted))
    }
    ##res$deviance <- -2 * res$logLik
    nvars <- NCOL(X) - 1 + length(unique(cluster))
    res$df.residual <- length(Y) - nvars
    res$n <- NROW(Y)
    res$aic <- res$deviance + 2 * nvars ## CHECK this !!
    res$boot <- TRUE
    res$call <- cl
    ##res$frail <- ifelse(res$frail > 999, Inf, res$frail)
    ##res$frail <- ifelse(res$frail < -999, -Inf, res$frail)
    if (!is.null(res$coefficients))
      names(res$coefficients) <- c(colnames(X))[-1]
    class(res) <- "glmmboot"
    res
}

