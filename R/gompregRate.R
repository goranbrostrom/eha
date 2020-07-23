gompregRate <- function(X, Y, strata, offset, init, control){
    ## Gompertz proportional hazards:
    ## Stratum i: h_i(t; beta) = a_i * exp(t * b_i) * exp(x*beta),
    ## or: h_i(t) = exp(shape_i + t * rate_i + x * beta)
    ## i = 1, ..., ns
    ## NOTE: Y is nn x 3!
    
    if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    nn <- NROW(X)
    ncov <- NCOL(X)
    if (NROW(Y) != nn) stop("Y NROW error")
    if (NCOL(Y) != 3) stop("Y NCOL error", )
    if (ncov){
        wts <- Y[, 2] - Y[, 1]
        means <- apply(X, 2, weighted.mean, w = wts)
        means <- apply(X, 2, mean) ## ??

        for (i in 1:ncov){
            X[, i] <- X[, i] - means[i]
        }

    }
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

    if (missing(init) || is.null(init)) # Should not happen except in direct call.
        init <- rep(0, ncov)
    if (length(init) != ncov) stop("Error in init")

    printlevel <- control$trace
    iter <- control$maxiter

    ## nstra <- c(0, cumsum(table(strata)))

    bdim <- ncov + 2 * ns

    dcumhaz_rate <- function(x, shape, rate){
        ## Gives the partial derivative of the cumulative hazards function
        ## wrt 'rate'. Deals also with the case 'rate = 0'.
        
        if (isTRUE(all.equal(rate, 0))){
            ret <- shape * x^2 / 2
        }else{
            eta <- exp(x * rate)
            ret <- (shape * rate * x * eta - shape * (eta - 1)) / rate^2
        }
        ret
    }

    dGomp <- function(beta){
        
        ## Calculates the first derivatives of a Gompertz regression (stratified)
        ## beta[1:ncov] = the beta coefficients.
        ## beta[ncov + 1], beta[ncov + 3], ... = rate[1, 2, ...] ("gamma")
        ## beta[ncov + 2], beta[ncov + 4], ... = shape[1, 2, ...] ("alpha")
        
        enter <- Y[, 1]
        exit <- Y[, 2]
        event <- Y[, 3]
        if (ncov){
            bz <- offset + X %*% beta[1:ncov]
        }else{
            bz <- offset
        }
        
        grad <- numeric(ncov + 2 * ns)
        
        for (i in 1:ns){
            T0 <- enter[strata == i]
            T <- exit[strata == i]
            D <- event[strata == i]
            z <- X[strata == i, ,drop = FALSE]
            
            ezb <- exp(bz[strata == i])
            shape <- exp(beta[ncov + 2 * i]) # beta[ncov+2i] = log(shape)
            rate <- beta[ncov + 2 * i - 1] 
            
            ##ezbH <- ezb * (exp(rate * T) - exp(rate * T0))
            ezbHr <- ezb * (dcumhaz_rate(T, shape, rate) - 
                               dcumhaz_rate(T0, shape, rate))
            
            ezbH <- ezb * 
                (Hgompertz(T, shape = shape, rate = rate, param = "rate") -
                     Hgompertz(T0, shape = shape, rate = rate, param = "rate"))
            ezbHs <- ezbH 
            
            ## shape:
            grad[ncov + 2 * i] <- sum(D) - sum(ezbHs)
            
            ##rate:
            grad[ncov + 2 * i - 1] <- sum(D * T) - sum(ezbHr)
            
            ## Covariates (if any):
            if (ncov){
                for (j in 1:ncov){
                    grad[j] <- grad[j] +
                        ##sum(D * z[, j]) - sum(ezbH * z[, j]) * exp(shape) / rate
                        sum(D * z[, j]) - sum(z[, j] * ezbH)
                }
            }
        }

        grad
    }    
    
    Fmin <- function(beta){
        
        total <- 0
        
        enter <- Y[, 1]
        exit <- Y[, 2]
        event <- Y[, 3]
        
        for (i in 1:ns){
            T0 <- enter[strata == i]
            T <- exit[strata == i]
            D <- event[strata == i]
            rate <- beta[ncov + 2 * i - 1] # Note; eg, log(scale) = rate!

            shape <- exp(beta[ncov + 2 * i])    # Note; beta = log(shape)!
            if (ncov){
                bz <- offset[strata == i] +
                    X[strata == i, , drop = FALSE] %*% beta[1:ncov]
                if (printlevel) cat("beta = ", beta[1:ncov], "\n")
            }else{
                bz <- offset[strata == i]
            }
            ebz <- exp(bz)
            H1 <- Hgompertz(T, shape = shape, rate = rate, param = "rate")
            H0 <- Hgompertz(T0, shape = shape, rate = rate, param = "rate")
            ## log(h):
            h <- hgompertz(T, shape = shape, rate = rate, 
                           param = "rate", log = TRUE)
            total <- total + sum(D * (h + bz)) - sum(ebz * (H1 - H0))
        }
        if (printlevel){
            cat("[Fmin]: total = ", total, "\n")
        }
        return(total)
    }

    
    ## First, fit the 'null' model:
    beta0 <- numeric(2 * ns)
    ncov.save <- ncov
    ncov <- 0
    
    ## Start values for NULL model:

    for (i in 1:ns){
        enter <- Y[strata == i, 1]
        exit <- Y[strata == i, 2]
        event <- Y[strata == i, 3]
        score <- offset[strata == i]## + X[strata == i, , drop = FALSE] %*% init
        beta0[(2 * i - 1):(2 * i)] <-
            gompstartRate(enter, exit, event, score, simple = FALSE) 
        if (printlevel){
            cat("[gompregRate null]: i = ", i, "start values = ", beta0[(2 * i - 1):(2 * i)], "\n")
        }
    }

    ## NULL model:
    res0 <- optim(beta0, Fmin, gr = dGomp,
                 method = "BFGS",
                 control = list(fnscale = -1, reltol = 1e-10, 
                                trace = printlevel),
                 hessian = FALSE)


    if (printlevel){
        cat("counts: ", res0$counts, "\n")
        cat("[res0] fprim = ", dGomp(res0$par), "\n")
    }
    
    ## Done; now the real thing:    
    
    ncov <- ncov.save

    beta <- numeric(bdim)
    
    if (ncov)
        beta[1:ncov] <- init  # Start values

    beta[(ncov + 1):length(beta)] <- res0$par # Ditto
    if (FALSE){
        ##Start values: # Already given above
        for (i in 1:ns){
            enter <- Y[strata == i, 1]
            exit <- Y[strata == i, 2]
            event <- Y[strata == i, 3]
            score <- offset[strata == i] + X[strata == i, , drop = FALSE] %*% init
            beta[(ncov + 2 * i - 1):(ncov + 2 * i)] <-
                gompstartRate(enter, exit, event, score, simple = FALSE) 
            if (printlevel){
                cat("[gompregRate]: i = ", i, "start values = ", 
                    beta, "\n")
                ##beta[(ncov + 2 * i - 1):(ncov + 2 * i)], "\n")
            }
        }
    }

    ##ncov.save <- ncov
    ##ncov <- 0
    ##l0 <- Fmin(beta[(ncov.save + 1):bdim])
    ##ncov <- ncov.save
    res <- optim(beta, Fmin, gr = dGomp,
                 method = "BFGS",
                 control = list(fnscale = -1, reltol = 1e-10,
                                trace = printlevel),
                 hessian = TRUE)
    if (res$convergence != 0) stop("[gompreg]: No convergence")
    coefficients <- res$par
    coef.names <- colnames(X)
    if (ns > 1){
        for (i in 1:ns){
            coef.names <- c(coef.names,
                            paste("rate", as.character(i), sep =":"),
                            paste("log(level)", as.character(i), sep =":"))
        }
        
    }else{
        coef.names <- c(coef.names,
                        "rate", "log(level)")
    }

    names(coefficients) <- coef.names
    fit <- list(coefficients = coefficients,
                loglik = c(res0$value, res$value) #c(res0$value, res$value)
                )
    fit$gradient <- dGomp(coefficients)
    fit$pfixed <- FALSE
    fit$var <- tryCatch(solve(-res$hessian), error = function(e) e)
    if (is.matrix(fit$var)){
        colnames(fit$var) <- coef.names
        rownames(fit$var) <- coef.names
    }
    fit$hessian <- res$hessian
    if (is.matrix(fit$hessian)){
        colnames(fit$hessian) <- coef.names
        rownames(fit$hessian) <- coef.names
    }

    ## Fixing the shape parameter(s) due to centering:
    ##cat(" Before: fit$coefficients = ", fit$coefficients, "\n")
    ##cat("means = ", means, "\n")
    ##cat("coefficients = ", fit$coefficients[1], "\n")
    if (ncov){
        dxy <- diag(2 * ns + ncov)
        for (i in 1:ns){
            row <- ncov + 2 * i
            shape.corr <- sum(means * fit$coefficients[1:ncov]) #/
            fit$coefficients[row] <- fit$coefficients[row] - shape.corr
            dxy[row, 1:ncov] <- -means
        }
        if (is.numeric(fit$var)) fit$var <- dxy %*% fit$var %*% t(dxy)
        fit$dxy <- dxy
    }
    if (is.matrix(fit$var)){
        colnames(fit$var) <- coef.names
        rownames(fit$var) <- coef.names
    }
    ##cat(" After: fit$coefficients = ", fit$coefficients, "\n")
    
    fit$n.strata <- ns
    fit$df <- ncov
    fit$fail <- FALSE # Optimist!
    fit$param <- "rate"
    fit
}
                      
