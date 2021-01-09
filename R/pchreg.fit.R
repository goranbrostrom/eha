pchreg.fit <- function(X, Y, cuts, offset, strata, init, control){
    ## Test of stratification of the pch-model ("piecewise constant")
    ## strata \in {1, ..., ns}
    
    if (NCOL(Y) == 2) Y <- cbind(rep(0, NROW(Y)), Y)
    
    splitSurv <- function(Y, cuts){
        ## New way of cutting: Time is limited by min(cuts), max(cuts)
        ## That is, pieces outside that interval are cut off.
        ## Contrary to the eha::SurvSplit and survival::survSplit fcns!
        
        n <- length(cuts) - 1 # No. of intervals
        if (n < 1) stop("Number of cuts must be at least 2.")
        indat <- cbind(Y, 1:NROW(Y), rep(-1, NROW(Y)))
        colnames(indat) <- c("enter", "exit", "event", "idx", "ivl")
        cuts <- sort(cuts)
        indat <- as.data.frame(indat)
        out <- vector(mode = "list", length = n)
        for (i in 1:n){
            out[[i]] <- age.window(indat, cuts[i:(i+1)])
            out[[i]]$ivl <- i
        }
        Y <- do.call(rbind, out)
        colnames(Y) <- colnames(indat)
        list(Y = Y[, 1:3],
             ivl = Y[, 5],
             idx = Y[, 4]
        )
    }
    
    nn <- NROW(Y)
    if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    if (NROW(X) != nn) stop("[pchreg.fit]: error in X")
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
    
    if (missing(cuts)){
        cuts <- c(min(Y[, 1]), max(Y[, 2]))
    }
    
    n.ivl <- length(cuts) - 1
    if (n.ivl <= 0) stop("Need at least one interval (two or zero cut points).")
    ## Some sane checks .....
    ##cuts <- c(0, cuts, Inf)
    if (control$trace){
        cat("[pchreg.fit]    ns = ", ns, "\n")
        cat("[pchreg.fit] n.ivl = ", n.ivl, "\n")
    }
    if (length(cuts)){
        split <- splitSurv(Y, cuts)
        
        strat <- strata[split$idx]
        offset <- offset[split$idx]
        X <- X[split$idx, ,drop = FALSE]
        T <- split$Y$exit - split$Y$enter
        d <- split$Y$event
        ivl <- split$ivl
    }else{
        strat <- strata
        T <- Y[, 2] - Y[, 1]
        d <- Y[, 3]
        ivl <- rep(1, NROW(Y))
    }
    
    Dtot <- sum(d)


    loglik0 <- function(){
        alpha <- matrix(0, nrow = ns, ncol = n.ivl)
        res <- -Dtot
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strat == j) & (ivl == i)
                D <- sum(d[indx])
                if (D > 0){
                    sumT <- sum(T[indx])
                    alpha[j, i] <- D / sumT 
                    res <- res + D * (log(D) - log(sumT))
                }
            }
        }
        ##cat("res = ", res, "\n")
        list(loglik = res, hazards = alpha)
    }
    
    loglik <- function(beta){
        ##cat("beta = ", beta, "\n") 
        zb <- offset + X %*% beta
        ##cat("zb = ", zb, "\n")
        tezb <- T * exp(zb)
        ##cat("tezb = ", tezb, "\n")
        ##res <- sum(d * zb) - Dtot
        res <- drop(d %*% zb) - Dtot
        ##cat("res.first = ", res, "\n")
        ##cat("res1 = ", res, "\n")
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strat == j) & (ivl == i)
                D <- sum(d[indx])
                ##cat("i = ", i, " j = ", j, " D = ", D, "\n")
                if (D > 0){
                    res <- res + D * (log(D) - log(sum(tezb[indx])))
                }
            }
        }
        ##cat("res = ", res, "\n")
        res
    }

    dloglik <- function(beta){
        zb <- offset + X %*% beta
        tezb <- T * exp(zb)

        res <- drop(d %*% X) # 10 feb 2014
        ##res <- numeric(ncov)
        ##for (j in seq_len(ncov)){
          ##  res[j] <- sum(d * X[, j]) 
        ##}
        ##cat("res1[1] = ", res[1], "\n")
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strat == j) & (ivl == i)
                D <- sum(d[indx])
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

    getAlpha <- function(beta){
        ## Note: f = (alpha * exp(zb)^d * exp(-t * alpha * exp(zb))
        tezb <- T * exp(X %*% beta)
        alpha <- matrix(0, nrow = ns, ncol = n.ivl)
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strat == j) & (ivl == i)
                D <- sum(d[indx])
                alpha[j, i] <- D / sum(tezb[indx])
            }
        }
    alpha
    } # end getAlpha

    fit <- loglik0()
    
    if (ncov){
        means <- colMeans(X)
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
    }else{# No covariates
        fit$loglik <- rep(fit$loglik, 2)
        fit$coefficients <- NULL
    }
    
    fit$df <- ncov
    fit$cuts <- cuts
    fit$fail <- FALSE # Optimism...
    cn <- character(n.ivl)
    if (n.ivl > 1){
        ##cn[1] <- paste("(.., ", cuts[1], "]", sep = "")
        ##cn[n.ivl] <- paste("(", cuts[n.ivl - 1], ", ...]", sep = "")
        ##if (n.ivl > 2){
            for (i in 1:n.ivl){
                cn[i] <- paste("(", cuts[i], ", ", cuts[i + 1], "]", sep = "")
            }
        ##}
        colnames(fit$hazards) <- cn
    }
    
    if (ns > 1) rownames(fit$hazards) <- name.s
    fit$n.strata <- ns

    fit
}
