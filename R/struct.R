#' @export
struct <- function(m, covars, dropx, X, exposure){
    isI <- logical(NCOL(X)) # Added Jan 2014; 2.4-0.
    isF <- logical(length(covars))
    isO <- logical(length(covars))
    if (length(covars)){
        for (i in 1:length(covars)){
            if (length(dropx)){
                if (is.logical(m[, -(dropx + 1)][, (i + 1)])){
                    m[, -(dropx + 1)][, (i + 1)] <-
                        as.factor(m[, -(dropx + 1)][, (i + 1)])
                }
                isF[i] <- is.factor(m[, -(dropx + 1)][, (i + 1)])## ||
                isO[i] <- is.ordered(m[, -(dropx + 1)][, (i + 1)])## ||
                ##is.logical(m[, -(dropx + 1)][, (i + 1)]) )
            }else{
                if (is.logical(m[, (i + 1)])){
                    m[, (i + 1)] <- as.factor(m[, (i + 1)])
                }
                isF[i] <- is.factor(m[, (i + 1)]) ##||
                isO[i] <- is.ordered(m[, (i + 1)]) ##||
                ## is.logical(m[, (i + 1)]) )
            }
        }
    }
    
    if (any(isF)){
        levels <- list()
        index <- 0
        for ( i in 1:length(covars) ){
            if (isF[i]){
                index <- index + 1
                if (length(dropx)){
                    ll <- levels(m[, -(dropx + 1)][, (i + 1)])
                    if (isO[i]){
                        levels[[i]] <- paste("order", seq_along(ll) - 1)
                    }else{
                        levels[[i]] <- ll
                    }
                }else{
                    ll <- levels(m[, (i + 1)])
                    if (isO[i]){
                        levels[[i]] <- paste("order", seq_along(ll) - 1)
                    }else{
                        levels[[i]] <- ll
                    }
                }
            }else{
                levels[[i]] <- NULL
            }
        }
    }else{
        levels <- NULL
    }

    if (any(isF)){ ## New; get isI: (Jan 2014; 2.4-0)
        indx <- 0
        for (i in 1:length(covars)){
            indx <- indx + 1
            if (isF[i]){
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
    ###
    ttr <- sum(exposure)
    w.means <- list()
    if (length(covars)){
        for (i in 1:length(covars)){
            nam <- covars[i]
            ##cat("nam = ", nam, "\n")
            col.m <- which(nam == names(m))
            ##cat("col.m = ", col.m, "\n")
            if (isF[i]){
                n.lev <- length(levels[[i]])
                w.means[[i]] <- numeric(n.lev)
                for (j in 1:n.lev){
                    who <- m[, col.m] == levels[[i]][j]
                    ##cat("sum(exposure[who]) = ", sum(exposure[who]), "\n")
                    w.means[[i]][j] <-
                      sum(exposure[who] ) / ttr ## * 100, if in per cent
                }
            }else{
                w.means[[i]] <- sum(exposure * m[, col.m]) / ttr
            }
        }
    }


    list(isI = isI, isF = isF, isO = isO, levels = levels,
         w.means = w.means, ttr = ttr)
}   
