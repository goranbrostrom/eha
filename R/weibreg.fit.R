#' Weibull regression
#' 
#' This function is called by \code{\link{weibreg}}, but it can also be
#' directly called by a user.
#' 
#' See \code{\link{weibreg}} for more detail.
#' 
#' @param X The design (covariate) matrix.
#' @param Y A survival object, the response.
#' @param strata A stratum variable.
#' @param offset Offset.
#' @param init Initial regression parameter values.
#' @param shape If positive, a fixed value of the shape parameter in the
#' Weibull distribution. Otherwise, the shape is estimated.
#' @param control Controls convergence and output.
#' @param center Should covariates be centered?
#' @return \item{coefficients}{Estimated regression coefficients plus estimated
#' scale and shape coefficients, sorted by strata, if present.} \item{var}{}
#' \item{loglik}{Vector of length 2. The first component is the maximized
#' loglihood with only scale and shape in the model, the second the final
#' maximum.} \item{score}{Score test statistic at initial values}
#' \item{linear.predictors}{Linear predictors for each interval.}
#' \item{means}{Means of the covariates} \item{conver}{TRUE if convergence}
#' \item{fail}{TRUE if failure} \item{iter}{Number of Newton-Raphson iterates.}
#' \item{n.strata}{The number of strata in the data.}
#' @author Göran Broström
#' @seealso \code{\link{weibreg}}
#' @keywords survival regression
#' @export weibreg.fit
weibreg.fit <- function(X, Y,
                        strata, offset,
                        init, shape,
                        control, center = TRUE){

    ##X <- scale(X, center = TRUE, scale = FALSE)
    nn <- NROW(X)
    ncov <- NCOL(X)

    if (ncov) means <- colMeans(X)
    
    if (missing(strata) || is.null(strata)){ 
        strata <- rep(1, nn)
        ns <- 1
    }else{
        strata <- as.integer(factor(strata))
        ns <- max(strata)
    }
    
    if (length(strata) != nn) stop("Error in stratum variable")
    if (missing(offset) || is.null(offset)) 
      offset <- rep(0, nn)
    
    if (missing(init) || is.null(init)) 
      init <- rep(0, ncov)
    if (length(init) != ncov) stop("Error in init")
    
    printlevel <- control$trace
    iter <- control$maxiter

    if (shape <= 0){
        
        bdim <- ncov + 2 * ns
        if (ns > 0){
            ord <- order(strata)
            X <- X[ord, , drop = FALSE]
            Y <- Y[ord, , drop = FALSE]
            offset <- offset[ord]
            nstra <- c(0, cumsum(table(strata)))
        }
        ##if (center) X <- scale(X, center = TRUE, scale = FALSE)
        fit <- .C("weibsup",
                  iter = as.integer(iter), #maxit on ip, actual on op.
                  as.double(control$eps),
                  as.integer(printlevel),
                                        #
                  as.integer(ns), # No. of strata
                  as.integer(nstra),
                  as.integer(nn),
                  as.integer(ncov),
                  as.integer(bdim),
                                        #
                  as.double(Y[, 1]),  ## 'enter'
                  as.double(Y[, 2]),  ## 'exit'
                  as.integer(Y[, 3]), ## 'event'
                                        #
                  as.double(t(scale(X, center = TRUE, scale = FALSE))),
                  ## NOTE transpose!
                  as.double(offset),
                                        #
                  as.double(init),     # 'start.beta'
                  beta = double(bdim), # results -->
                  lambda = double(ns),
                  lambda.sd = double(ns),
                  shape = double(ns),  ## "p"
                  shape.sd = double(ns),
                                        #
                  loglik = double(2), # [1] == start, [2] == maximized
                  dloglik = double(bdim),
                  variance = double(bdim * bdim),
                  sctest = double(1),
                                        #
                  conver = integer(1),
                  fail = integer(1),
                  ## DUP = FALSE,
                  PACKAGE = "eha"
                  )
        
        if (fit$fail) return(list(fail = fit$fail,
                                  n.strata = ns,
                                  value = fit$beta[fit$fail])
                             )
        
        dxy <- diag(2 * ns + ncov)
        for (i in 1:ns){ ## Really a HACK ??!!!!!!!!!!!!!!!
            row <- ncov + 2 * i - 1
            col <- row + 1
            fit$beta[row] <- -fit$beta[row]
            if (ncov){
                pi.hat <- exp(fit$beta[col])
                scale.corr <- sum(means * fit$beta[1:ncov]) /
                  pi.hat
                fit$beta[row] <- fit$beta[row] + scale.corr
                dxy[row, 1:ncov] <- means / pi.hat
                dxy[row, col] <- -scale.corr
            }
            dxy[row, row] <- -1
        }
        coef.names <- colnames(X)
        if (ns > 1){
            for (i in 1:ns){
                coef.names <- c(coef.names,
                                paste("log(scale)", as.character(i), sep =":"),
                                paste("log(shape)", as.character(i), sep =":"))
            }
            
        }else{
            coef.names <- c(coef.names,
                            "log(scale)", "log(shape)")
        }
            
        
        fit$shape.fixed <- FALSE
    }else{  ## "Exponential" regression:
        if (ns >= 2) warning("'strata' is not useful for 'fixed shape' regression.\n Include stratum variable as a factor in the model instead.")
        bdim <- ncov + 1
        ##if (center) X <- scale(X, center = TRUE, scale = FALSE)
        
        fit <- .C("expsup",
                  iter = as.integer(iter), #maxit on ip, actual on op.
                  as.double(control$eps),
                  as.integer(printlevel),
                                        #
                  as.integer(nn),
                  as.integer(ncov),
                  as.integer(bdim),
                                        #
                  as.double(Y[, 1]),  ## 'enter'
                  as.double(Y[, 2]),  ## 'exit'
                  as.integer(Y[, 3]), ## 'event'
                                        #
                  as.double(t(scale(X, center = TRUE, scale = FALSE))),
                  as.double(offset),
                  as.double(shape), ## "p"
                                        #
                  as.double(init),     # 'start.beta'
                  beta = double(bdim), # results -->
                  lambda = double(1),
                  lambda.sd = double(1),
                                        #
                  loglik = double(2), # [1] == start, [2] == maximized
                  dloglik = double(bdim),
                  variance = double(bdim * bdim),
                  sctest = double(1),
                                        #
                  conver = integer(1),
                  fail = integer(1),
                  ## DUP = FALSE,
                  PACKAGE = "eha"
                  )
        if (fit$fail) return(list(fail = fit$fail,
                                  n.strata = ns,
                                  value = fit$beta[fit$fail])
                             )
        fit$shape.fixed <- TRUE
        fit$shape <- shape
        fit$shape.sd <- NULL  ## Not necessary!?!?
        fit$beta[bdim] <- -fit$beta[bdim] ## To get "1 / lambda"!
        dxy <- diag(bdim)
        if (ncov){
            dxy[bdim, 1:ncov] <- means / shape
            scale.corr <- sum(means * fit$beta[1:ncov]) / shape
            fit$beta[bdim] <-
              fit$beta[bdim] + scale.corr
        }
        dxy[bdim, bdim] <- -1
        coef.names <- c(colnames(X), "log(scale)")
        
        ## Note; this is really a "hack"!!!!!!!!!!!!!!!
    }
    
    ##cat("done!\n")
    
    if (!fit$fail){
        var <- dxy %*% matrix(fit$variance, bdim, bdim) %*% t(dxy)
        colnames(var) <- rownames(var) <- coef.names
    }
    else
      var <- NULL

    coefficients <- fit$beta
    names(coefficients) <- coef.names 
    
    list(coefficients = coefficients,
         df = ncov,
         var = var,
         loglik = fit$loglik,
         score = fit$sctest,
         conver = fit$conver,
         fail = fit$fail,
         iter = fit$iter,
         n.strata = ns,
         shape = fit$shape
         )
    
}
