tpchreg.fit <- function(X, count, exposure, offset, weights, strata, time){
    ## NOTE: 'time' must be a factor here.
    
    if (!is.factor(time)){
        cat("time:", time[1:5], "\n")
        stop("[tpchreg.fit] time must be a factor here.")
    } 
    ##print(str(time))
    if (!missing(weights)){
        count <- count * weights
        exposure <- exposure * weights   
    }
    
    ivl <- as.integer(time) ## time must be a factor here!
    ##cat("unique(time) = ", unique(time), "\n")
    nn <- length(count)
        if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    if (NROW(X) != nn) stop("[pchreg2]: error in X")
    if (missing(strata) || is.null(strata)){
        strata <- rep(1, nn)
        ns <- 1
    }else{
        strata <- factor(strata)
        name.s <- levels(strata)
        strata <- as.integer(strata)
        ns <- max(strata)
    }
    ns <- length(unique(strata))
    ncov <- NCOL(X) # Assume here that ncov >= 0!
    init <- rep(0, ncov) ## FIX THIS hack!!
    n.ivl <- length(unique(time))
    Dtot <- sum(count) 
#########
    loglik0 <- function(){
        alpha <- matrix(0, nrow = ns, ncol = n.ivl)
        res <- -Dtot
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                ##cat(" i = ", i, "j = ", j, "\n")
                indx <- (strata == j) & (ivl == i)
                ##cat("indx = ", indx, "\n")
                ##cat("ivl = ", ivl, "\n")
                ##cat("sum(indx) = ", sum(indx), "\n")
                D <- sum(count[indx])
                ##cat("D = ", D, "\n")
                if (D > 0){
                    sumT <- sum(exposure[indx])
                    alpha[j, i] <- D / sumT 
                    res <- res + D * (log(D) - log(sumT))
                }
            }
        }
        ##cat("res = ", res, "\n")
        ##

        cbind(strata, ivl, count)
        list(loglik = res, hazards = alpha)
    }
#########
   loglik <- function(beta){
        ##cat("beta = ", beta, "\n") 
        zb <- offset + X %*% beta
        ##cat("zb = ", zb, "\n")
        tezb <- exposure * exp(zb) 
        ##cat("tezb = ", tezb, "\n")
        ##res <- sum(d * zb) - Dtot
        res <- drop(count %*% zb) - Dtot
        ##cat("res.first = ", res, "\n")
        ##cat("res1 = ", res, "\n")
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strata == j) & (ivl == i)
                D <- sum(count[indx])
                ##cat("i = ", i, " j = ", j, " D = ", D, "\n")
                if (D > 0){
                    res <- res + D * (log(D) - log(sum(tezb[indx])))
                }
            }
        }
        ##cat("res = ", res, "\n")
        res
    }
############
   dloglik <- function(beta){
        zb <- offset + X %*% beta
        tezb <- exposure * exp(zb)

        res <- drop(count %*% X) # 10 feb 2014
        ##res <- numeric(ncov)
        ##for (j in seq_len(ncov)){
          ##  res[j] <- sum(d * X[, j]) 
        ##}
        ##cat("res1[1] = ", res[1], "\n")
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strata == j) & (ivl == i)
                D <- sum(count[indx])
                ##cat("i = ", i, " j = ", j, " D = ", D, "\n")
                if (D > 0){
                    Stezb <- sum(tezb[indx]) 
                    for (j in seq_len(ncov)){ # Do smarter!
                        res[j] <- res[j] -
                            D * sum(X[indx, j] * tezb[indx]) /
                             Stezb
                    }
                }
            } # end for (j in ...
        } # end for (i in ...
        ##cat("res[1] = ", res[1], "\n")
        res
    } # end dloglik
############
    getAlpha <- function(beta){
        ## Note: f = (alpha * exp(zb)^d * exp(-t * alpha * exp(zb))
        tezb <- exposure * exp(X %*% beta)
        alpha <- matrix(0, nrow = ns, ncol = n.ivl)
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strata == j) & (ivl == i)
                D <- sum(count[indx])
                if (D > 0){
                    alpha[j, i] <- D / sum(tezb[indx])
                }
            }
        }
    alpha
    } # end getAlpha
############
     fit <- loglik0()
    
    if (ncov){
        ##means <- colMeans(X)
        means <- drop((exposure / sum(exposure)) %*% X)
        for (i in seq_len(ncov)){
            X[, i] <- X[, i] - means[i]
        }
        beta <- init
        res <- optim(beta, loglik, dloglik,
                     method = "BFGS", hessian = TRUE, 
                     control = list(fnscale = -1, reltol = 1e-10))
        beta <- res$par
        fit$gradient <- dloglik(beta)
        ##cat("score = ", deriv, " at solution\n")
        fit$loglik <- c(fit$loglik, res$value)
        fit$coefficients <- beta
        names(fit$coefficients) <- colnames(X)
        fit$var <- solve(-res$hessian)
        fit$hazards <- getAlpha(beta)
        fit$hazards <- fit$hazards * exp(-drop(means %*% fit$coefficients))
        colnames(fit$var) <- rownames(fit$var) <- colnames(X)
        fit$nullModel <- FALSE
        ##fit$w.means <- means ## MUST be fixed: This is WRONG!!
    }else{# No covariates
        fit$loglik <- rep(fit$loglik, 2)
        fit$coefficients <- NULL
        fit$nullModel <- TRUE
    }
############
    fit$df <- ncov
    ##fit$cuts <- cuts
    fit$fail <- FALSE # Optimism...
    ## cn <- character(n.ivl)
    ## if (n.ivl > 1){
    ##     cn[1] <- paste("(.., ", cuts[1], "]", sep = "")
    ##     cn[n.ivl] <- paste("(", cuts[n.ivl - 1], ", ...]", sep = "")
    ##     if (n.ivl > 2){
    ##         for (i in 2:(n.ivl - 1)){
    ##             cn[i] <- paste("(", cuts[i-1], ", ", cuts[i], "]", sep = "")
    ##         }
    ##     }
    ##     colnames(fit$hazards) <- cn
    ## }
    
    if (ns > 1) rownames(fit$hazards) <- name.s
    fit$n.strata <- ns

    fit

}
