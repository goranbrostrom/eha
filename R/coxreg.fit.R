#' Cox regression
#' 
#' Called by \code{\link{coxreg}}, but a user can call it directly.
#' 
#' \code{rs} is dangerous to use when NA's are present.
#' 
#' @param X The design matrix.
#' @param Y The survival object.
#' @param rs The risk set composition. If absent, calculated.
#' @param weights Case weights; time-fixed or time-varying.
#' @param t.offset Case offset; time-varying.
#' @param strats The stratum variable. Can be absent.
#' @param offset Offset. Can be absent.
#' @param init Start values. If absent, equal to zero.
#' @param max.survs Sampling of risk sets? If so, gives the maximum number of
#' survivors in each risk set.
#' @param method Either "efron" (default) or "breslow".
#'
#' @param boot Number of bootstrap replicates. Defaults to FALSE, no
#' bootstrapping.
#' @param efrac Upper limit of fraction failures in 'mppl'.
#' @param calc.martres Should martingale residuals be calculated?
#' @param control See \code{\link{coxreg}}
#' @param verbose Should Warnings about convergence be printed?
#' @param calc.hazards Deprecated. See \code{\link{coxreg}}.
#' @param center Deprecated. See \code{\link{coxreg}}.
#'  
#' @return A list with components \item{coefficients}{Estimated regression
#' parameters.} \item{var}{Covariance matrix of estimated coefficients.}
#' \item{loglik}{First component is value at \code{init}, second at maximum.}
#' \item{score}{Score test statistic, at initial value.}
#' \item{linear.predictors}{Linear predictors.} \item{residuals}{Martingale
#' residuals.} \item{hazard}{Estimated baseline hazard. At value zero of
#' 'design' variables.} \item{means}{Means of the columns of the design
#' matrix.} \item{bootstrap}{The bootstrap replicates, if requested on input.}
#' \item{conver}{\code{TRUE} if convergence.} \item{f.conver}{TRUE if variables
#' converged.} \item{fail}{\code{TRUE} if failure.} \item{iter}{Number of
#' performed iterations.}
#' @note It is the user's responsibility to check that indata is sane.
#' @author Göran Broström
#' @seealso \code{\link{coxreg}}, \code{\link{risksets}}
#' @keywords survival regression
#' @examples
#' 
#'  X <- as.matrix(data.frame(
#'                 x=     c(0, 2,1,4,1,0,3),
#'                 sex=   c(1, 0,0,0,1,1,1)))
#'  time <- c(1,2,3,4,5,6,7)
#'  status <- c(1,1,1,0,1,1,0)
#'  stratum <- rep(1, length(time))
#' 
#'  coxreg.fit(X, Surv(time, status), strats = stratum, max.survs = 6,
#'      control = list(eps=1.e-4, maxiter = 10, trace = FALSE))
#' 
#' @export
coxreg.fit <- function(X, Y, rs, weights, t.offset = NULL,
                       strats, offset, init, max.survs,
                       method = "efron",
                       boot = FALSE, efrac = 0,
                       calc.martres = TRUE,
                       control, verbose = TRUE, calc.hazards = NULL,
                       center = NULL){

    if (!missing(center)){
      warning("argument 'center' is deprecated.")
    }
    if (!missing(calc.hazards)){
      warning("argument 'calc.hazards' is deprecated")
    }
  
    nn <- NROW(Y)
    if (is.matrix(X)){
        ncov <- NCOL(X)
    }else{
        if (length(X)) NCOL <- 1
        else NCOL <- 0
    }

    if (missing(strats) || is.null(strats))
      strats <- rep(1, nn)

    if (missing(rs) || is.null(rs)){
        rs <- risksets(Y, strata = strats, max.survs)
    }

    if (missing(weights) || is.null(weights)){
        weights <- rep(1, length(rs$riskset))
    }else{
        if (length(weights) == nn){
            weights <- weights[rs$riskset]
        }else if (length(weights) != length(rs$riskset)){
            stop("weights: length error")
        }
    }

    if (missing(t.offset) || is.null(t.offset)){
        t.offset <- rep(0, length(rs$riskset))
    }else{
        if (length(t.offset) != length(rs$riskset)){
            cat("length(t.offset) = ", length(t.offset), "\n")
            cat("length(rs$riskset) = ", length(rs$riskset), "\n")
            stop("t.offset: length error")
        }
    }

    if (max(rs$riskset) > nn) stop("Riskset does not match data")

    if (missing(offset) || is.null(offset)){
        offset <- t.offset
    }else{
        offset <- t.offset + offset[rs$riskset]
    }

    if (missing(init) || is.null(init))
        init <- rep(0, ncov)

    if (missing(control)){
        control <- list(eps=1.e-8, maxiter = 10, trace = FALSE)
    }else{
        if (!is.numeric(control$eps)){
            stop("Error in control = list(eps = ...) ")
        }else{
            if (control$eps <= 0) stop("control$eps must be strictly positive")
        }
        if (!is.numeric(control$maxiter)){
            stop("Error in control = list(maxiter = ...) ")
        }else{
            if (control$maxiter < 0) stop("control$maxiter must be positive")
        }
        if (!is.logical(control$trace)) stop("control$trace must be logical")
    }

    nullModel <- (ncov == 0)

    if (nullModel){
        ## faking a simple model with no iterations
        ncov <- 0
        X <- matrix(0, ncol = 1, nrow = nn)
        control$maxiter <- 0
        init <- 0
        means <- NULL
    }else{
        means = apply(X, 2, mean)
    }

    printlevel <- control$trace
    ## NOTE: silent == TRUE ===> printlevel = 0
    iter <- control$maxiter
    if (method[1] == "efron")
      meth <- 0
    else if (method[1] == "breslow")
      meth <- 1
    else if (method[1] == "mppl")
      meth <- 2
    else if (method[1] == "ml")
      meth <- 3
    else
      stop(paste("Unknown method", as.character(method[1])))

    boot <- abs(as.integer(boot))

    ##if (center) X <- scale(X, center = TRUE, scale = FALSE)
    X <- scale(X, center = TRUE, scale = FALSE) # Always center!
    ## if (!nullModel){
        fit <- .C("sup",
                  as.integer(meth),
                  iter = as.integer(iter), #maxit on input, actual on output
                  as.double(control$eps),
                  as.integer(printlevel),
                                        #
                  as.integer(sum(rs$n.events)), ## total No. of events
                  as.integer(sum(rs$antrs)),  ## total No. of risksets
                  as.integer(length(rs$antrs)), # No. of strata
                                        #
                  as.integer(rs$antrs),
                  as.integer(rs$n.events),
                  as.integer(rs$size),
                  as.double(weights),
                                        #
                  as.integer(length(rs$riskset)), # Sum of risk set sizes.
                  as.integer(rs$eventset),
                  as.integer(rs$riskset),
                                        #
                  as.integer(nn),
                  as.integer(ncov),
                  ## Note here; X is transposed! From 2007-04-16 (0.99)
                  as.double(t(X)),
                  as.double(offset),
                                        #
                  as.double(init),     # 'start.beta'
                  boot = as.integer(boot),
                  as.double(efrac),
                  beta = double(ncov * (1 + boot)),
                  sd.beta = double(ncov * (1 + boot)),
                                        #
                  loglik = double(2), # [1] == start, [2] == maximized
                  variance = double(ncov * ncov),
                  sctest = double(1),
                  ##
                  hazard = double(sum(rs$antrs)),
                  conver = integer(1),
                  f.conver = integer(1),
                  fail = integer(1),
                  ## DUP = FALSE,
                  PACKAGE = "eha")

    if (FALSE){ ## NO!!!! 20110105; # YES!!! 20110103. Not for the moment...
        score.means <- exp(sum(means * fit$beta[1:ncov]))
        haz.mean <- 1 - (1 - fit$hazard)^score.means
    }else{
        haz.mean <- fit$hazard
    }
    hazards <- list()
    stopp <- cumsum(rs$antrs)
    startt <- c(1, 1 + stopp[-length(rs$antrs)])
    for (i in 1:length(rs$antrs)){
        hazards[[i]] <- cbind(rs$risktimes[startt[i]:stopp[i]],
                              haz.mean[startt[i]:stopp[i]])
    }
    class(hazards) <- "hazdata"
    bootstrap <- NULL
    boot.sd <- NULL
    if (boot & (fit$fail == 0)){

        bootstrap <- matrix(fit$beta[(ncov + 1):((boot + 1) * ncov)],
                            nrow = ncov, ncol = boot)
        boot.sd <- matrix(fit$sd.beta[(ncov + 1):((boot + 1) * ncov)],
                          nrow = ncov, ncol = boot)
        fit$beta <- fit$beta[1:ncov]
        fit$sd.beta <- fit$sd.beta[1:ncov]
    }
    if (FALSE) { # if nullModel
        X <- matrix(0, ncol = 0, nrow = nn)
        fit <- list(beta = numeric(0),
                    conver = TRUE,
                    fail = FALSE)
        ncov <- 0
        bootstrap <- NULL
        boot.sd <- NULL

    }
    if (fit$fail){
        out <- paste("Singular hessian; suspicious variable No. ",
                     as.character(fit$fail), ":\n",
                     colnames(X)[fit$fail], sep = "")
        if (verbose) warning(out)## New
        return(fit)## 19 February 2007.
    }else if (!fit$conver){
        ##fit$conver <- 1 Removed 19 February 2007
        if (!fit$f.conver){
            if (verbose)
              warning("Did not converge")
        }else{
            if (verbose)
              warning("log likelihood converged, but not variables")
        }
    }

    if (FALSE){ ##(calc.hazards && (!fit$fail)){
        score <- exp(X %*% fit$beta)
        hazard <- .Fortran("hazards",
                           as.integer(sum(rs$antrs)),  ## total No. of risksets
                           as.integer(length(rs$antrs)), # No. of strata
                           ##
                           as.integer(rs$antrs),
                           as.integer(rs$n.events),
                           as.integer(rs$size),
                           ##
                           as.integer(length(rs$riskset)),
                           ## is Sum of risk set sizes.
                           as.integer(rs$riskset),
                           ##
                           as.integer(nn),
                           ##
                           as.double(score),
                           hazard = double(sum(rs$antrs)),
                           ##
                           ##DUP = FALSE,
                           PACKAGE = "eha"
                           )$hazard
        ##}else{ #if not calc.hazards or fail
        hazards <- NULL
    }

    if (calc.martres && (!nullModel)){
        score <- exp(X %*% fit$beta)
        resid <- .Fortran("martres",
                          as.integer(sum(rs$antrs)),
                          as.integer(length(rs$antrs)),
                          ##
                          as.integer(rs$antrs),
                          as.integer(rs$n.events),
                          as.integer(rs$size),
                          ##
                          as.integer(length(rs$riskset)), # Sum of risk set sizes.
                          as.integer(rs$riskset),
                          ##
                          as.integer(nn),
                          ##
                          as.double(score),  # 'score'
                          as.double(fit$hazard),
                          resid = double(nn),
                          ##DUP = FALSE,
                          PACKAGE = "eha"
                          )$resid
    }else{
        resid <- NULL
    }

    if ((!fit$fail) && (!nullModel))
      var <- matrix(fit$variance, ncov, ncov)
    else
      var <- NULL


    list(coefficients = fit$beta,
         df = length(fit$beta),
         n = NROW(X),
         sd = fit$sd.beta,
         var = var,
         loglik = fit$loglik,
         score = fit$sctest,
         linear.predictors = X %*% fit$beta,
         residuals = resid,
         noOfRisksets = rs$antrs,
         ##risktimes = rs$risktimes,
         ##hazard = hazard,
         hazards = hazards,
         means = means,
         bootstrap = bootstrap,
         boot.sd = boot.sd,
         conver = fit$conver,
         f.conver = fit$f.conver,
         fail = fit$fail,
         iter = fit$iter
         )
}
