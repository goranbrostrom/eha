#' Loglihood function (partial likelihood) of a Cox regression
#' 
#' Calculates minus the log likelihood function and its first and second order
#' derivatives for data from a Cox regression model. It is used by \code{coxreg}
#' if the argument \code{coxph = FALSE}
#' 
#' Note that the function returns log likelihood, score vector and minus
#' hessian, i.e. the observed information. The model is 
#' 
#' @param beta Regression parameters
#' @param X The design (covariate) matrix.
#' @param offset Offset.
#' @param rs Risk set created by \code{risksets(..., collate_sets = TRUE)}
#' @param what what = 0 means only loglihood, 1 means score vector as well, 2
#' loglihood, score and hessian.
#' @return A list with components 
#' \item{loglik}{The log likelihood.}
#' \item{dloglik}{The score vector. Nonzero if \code{what >= 1}}
#' \item{d2loglik}{The hessian. Nonzero if \code{ord >= 2}}
#' @author Göran Broström
#' @seealso \code{\link{coxreg}}
#' @keywords survival distribution
#' @export coxfunk
##
coxfunk <- function(beta, X, offset, rs, what = 2){
    ##  Partial likelihood calculation.
    ## X: a design matrix nxk
    ## beta vector of length k
    ## offset: vector of length n
    ## list with risk sets.
    ##
    ## Note: 'rs' contains stratum information.
    
    lin.pred <- X %*% beta + offset
    score <- exp(lin.pred)
    
    p <- NCOL(X)
    loglik <- 0
    dloglik <- numeric(p)
    d2loglik <- matrix(0, nrow = p, ncol = p)
    nn <- length(rs)

    for (i in 1:nn){
        who <- rs[[i]]$risk
        antevents <- length(rs[[i]]$events)
        size <- length(who)
        ## Hack:
        weights <- rep(1, size)
        ##offs <- offset[who]
        lin <- lin.pred[who]
        xt <- X[who, ]
        e_frac <- rs[[i]]$sfrac
        ##cat("Entering C: ")
        res <- .C("breslow_rs2",
                  as.integer(what),
                  as.integer(antevents),
                  as.integer(size),
                  as.double(weights),
                  as.double(t(xt)),
                  as.double(lin),
                  as.integer(p),
                  as.double(beta),
                  as.double(e_frac),
                  ## Return:
                  loglik = double(1),
                  dloglik = double(p),
                  d2loglik = double(p * p),
                  PACKAGE = "eha")
        ##cat(" dloglik = ", res$dloglik, "\n")
        loglik <- loglik + res$loglik
        dloglik <- dloglik + res$dloglik
        d2loglik <- d2loglik + res$d2loglik
    }

    list(loglik = loglik, dloglik = dloglik, d2loglik = d2loglik)
}
                   
