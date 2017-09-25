#' Generalized Linear Models with random intercept
#' 
#' Fits GLMs with random intercept by Maximum Likelihood and numerical
#' integration via Gauss-Hermite quadrature.
#' 
#' The integrals in the log likelihood function are evaluated by the Laplace
#' approximation (default) or Gauss-Hermite quadrature. The latter is now fully
#' adaptive; however, only approximate estimates of variances are available for
#' the Gauss-Hermite (n.points > 1) method.
#' 
#' For the binomial families, the response can be a two-column matrix, see the
#' help page for glm for details.
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
#' @param weights Case weights. Defaults to one.
#' @param cluster.weights Cluster weights. Defaults to one.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action See glm.
#' @param start.coef starting values for the parameters in the linear
#' predictor.  Defaults to zero.
#' @param start.sigma starting value for the mixing standard deviation.
#' Defaults to 0.5.
#' @param fix.sigma Should sigma be fixed at start.sigma?
#' @param x If TRUE, the design matrix is returned (as x).
#' @param offset this can be used to specify an a priori known component to be
#' included in the linear predictor during fitting.
#' @param prior Which "prior" distribution (for the random effects)? Possible
#' choices are "gaussian" (default), "logistic", and "cauchy". For the poisson
#' family, it is possible to use the conjugate "gamma" prior, which avoids
#' numerical integration.
#' @param control Controls the convergence criteria. See
#' \code{\link{glm.control}} for details.
#' @param method There are two choices "Laplace" (default) and "ghq"
#' (Gauss-Hermite).
#' @param n.points Number of points in the Gauss-Hermite quadrature. If
#' n.points == 1, the Gauss-Hermite is the same as Laplace approximation. If
#' \code{method} is set to "Laplace", this parameter is ignored.
#' @param boot Do you want a bootstrap estimate of cluster effect? The default
#' is \emph{No} (\code{boot = 0}). If you want to say yes, enter a positive
#' integer here. It should be equal to the number of bootstrap samples you want
#' to draw. A recomended absolute \emph{minimum value} is \code{boot = 2000}.
#' @return The return value is a list, an object of class 'glmmML'. The
#' components are: \item{boot}{No. of boot replicates}
#' \item{converged}{Logical} \item{coefficients}{Estimated regression
#' coefficients} \item{coef.sd}{Their standard errors} \item{sigma}{The
#' estimated random effects' standard deviation} \item{sigma.sd}{Its standard
#' error} \item{variance}{The estimated variance-covariance matrix. The last
#' column/row corresponds to the standard deviation of the random effects
#' (\code{sigma})} \item{aic}{AIC} \item{bootP}{Bootstrap p value from testing
#' the null hypothesis of no random effect (sigma = 0)}
#' \item{deviance}{Deviance} \item{mixed}{Logical} \item{df.residual}{Degrees
#' of freedom} \item{cluster.null.deviance}{Deviance from a glm with no
#' clustering. Subtracting \code{deviance} gives a test statistic for the null
#' hypothesis of no clustering. Its asymptotic distribution is a symmetric
#' mixture a constant at zero and a chi-squared distribution with one df. The
#' printed p-value is based on this.} \item{cluster.null.df}{Its degrees of
#' freedom} \item{posterior.modes}{Estimated posterior modes of the random
#' effects} \item{terms}{The terms object} \item{info}{From hessian inversion.
#' Should be 0. If not, no variances could be estimated. You could try fixing
#' sigma at the estimated value and rerun.} \item{prior}{Which prior was used?}
#' \item{call}{The function call} \item{x}{The design matrix if asked for,
#' otherwise not present}
#' @note The optimization may not converge with the default value of
#' \code{start.sigma}. In that case, try different start values for sigma. If
#' still no convergence, consider the possibility to fix the value of sigma at
#' several values and study the profile likelihood.
#' @author Göran Broström
#' @seealso \code{\link{glmmboot}}, \code{\link{glm}}, \code{\link{optim}},
#' \code{lmer} in the package \code{lme4} and \code{glmmPQL} in the package
#' \code{MASS}.
#' @references Broström (2003). Generalized linear models with random
#' intercepts. \url{http://www.stat.umu.se/forskning/reports/glmmML.pdf}
#' @keywords regression
#' @examples
#' 
#' id <- factor(rep(1:20, rep(5, 20)))
#' y <- rbinom(100, prob = rep(runif(20), rep(5, 20)), size = 1)
#' x <- rnorm(100)
#' dat <- data.frame(y = y, x = x, id = id)
#' glmmML(y ~ x, data = dat, cluster = id)
#' 
#' @export glmmML
glmmML <- function(formula,
                   family = binomial,
                   data,
                   cluster,
                   weights,
                   cluster.weights,
                   subset,
                   na.action,
                   offset,
                   prior = c("gaussian", "logistic", "cauchy", "gamma"),
                   start.coef = NULL,
                   start.sigma = NULL,
                   fix.sigma = FALSE,
                   x = FALSE, # Should design matrix be returned?
                   control = list(epsilon = 1.e-8,
                       maxit = 200, trace = FALSE),
                   method = c("Laplace", "ghq"),
                   n.points = 8,
                   boot = 0){

    method <- method[1]
    if (method == "laplace") method <- "Laplace"
    if (method == "GHQ") method <- "ghq"
    if (!(method %in% c("Laplace", "ghq"))) stop("Wrong method")
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

    ## don't use this for the moment!
    ## Not used:
    
    method <- as.numeric(method[1] == "Laplace")
    if (method) n.points <- 1

    if (is.character(family)) 
        family <- get(family)
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
    
    a.prior <- prior[1]
    if (!(a.prior %in% c("gaussian", "logistic", "cauchy"))){
        if ((a.prior == "gamma") && (family$family == "binomial")){
            stop("The gamma prior only works with Poisson responses")
        }else if (a.prior != "gamma"){
            stop("Prior distribution not known")
        }
    }

    if (a.prior == "gaussian") prior <- 0
    else if (a.prior == "logistic") prior <- 1
    else if (a.prior == "gamma") prior <- 3
    else prior <- 2 # Cauchy
    
    ## 'gaussian' is the default
    
    cl <- match.call()

    if (missing(data))
        data <- environment(formula)
    
    mf <- match.call(expand.dots = FALSE)
    ## get a copy of the call; result: a list.
    
    mf$family <- mf$start.coef <- mf$start.sigma <- mf$fix.sigma <- NULL
    mf$weights <- mf$cluster.weights <- NULL
    mf$control <- mf$maxit <- mf$boot <- NULL
    mf$n.points <- mf$method <- mf$prior <- NULL

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


    ##    return(clus)
    
    ##if (NCOL(Y) >  1) stop("Response must be univariate")
    
    if (!is.null(offset) && length(offset) != NROW(Y)) 
        stop(paste("Number of offsets is", length(offset), ", should equal", 
                   NROW(Y), "(number of observations)"))
    
    
    if (missing(weights)) weights <- rep.int(1, NROW(Y))
    if (any(weights < 0)) stop("negative weights not allowed")

    if (missing(cluster.weights))
      cluster.weights <- rep.int(1, length(cluster))
    if (any(cluster.weights < 0)) stop("negative cluster weights not allowed")

    if (n.points <= 0) n.points <- 1 # Should give 'Laplace'(?)
    fit <- glmmML.fit(X, Y,
                      weights,
                      cluster.weights,
                      start.coef,
                      start.sigma,
                      fix.sigma,
                      cluster,
                      offset,
                      family,
                      method,
                      n.points,
                      control,
                      intercept = ( attr(mt, "intercept") > 0),
                      boot,
                      prior # gaussian by default
                      )
    
    if (!fit$convergence)
      warning("'vmmin' did not converge. Increase 'maxit'?")
##    if (fit$info) return(list(info = fit$info,
##                              convergence = fit$convergence,
##                              sigma = fit$sigma,
##                              coefficients = fit$beta,
##                              deviance = fit$deviance)
##                         )
    bdim <- p + 1
    res <- list()
    res$boot <- boot
    res$converged <- as.logical(fit$convergence)
    res$coefficients <- fit$beta
    res$coef.sd <- fit$beta.sd
    res$sigma <- abs(fit$sigma) # Note 
    res$sigma.sd <- fit$sigma.sd
    ## For the time being: Show the attained max!!!!!!!!!!!!!
    ##if (fit$cluster.null.deviance <= fit$deviance){
    ##      res$sigma = 0
    ##      res$sigma.sd = NA
    ##  }
    res$variance <- fit$variance
    res$aic <- fit$aic
    names(res$coef.sd) <- names(res$coefficients)
    
    res$bootP <- fit$bootP
    res$deviance <- fit$deviance
    res$df.residual <- fit$df.residual
    res$cluster.null.deviance <- fit$cluster.null.deviance
    res$cluster.null.df <- fit$cluster.null.df
    res$posterior.modes <- fit$post.mode
##    res$posterior.means <- fit$post.mean
    res$prior <- a.prior
    res$terms <- mt
    res$info <- fit$info # From inverting the hessian! Should be zero.
    res$call <- cl
    if (x) res$x <- X
    names(res$coefficients) <- c(colnames(X))
    class(res) <- "glmmML"
    res
}

