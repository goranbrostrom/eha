coxfun0 <- function(Y, strata = NULL, method = "efron"){
    rs <- risksets(Y, strata = strata, members = FALSE)
    d <- rs$n.events
    n <- rs$size
    r <- length(d)
    if (method == "efron"){
        ll <- 0
        for (i in 1:r){
            tmp <- 0
            for (j in 1:d[i]){
                tmp <- tmp + log(n[i] - j + 1)
            }
            ll <- ll - tmp
        }
    }else{
        ll <- -sum(d * log(n))
    }
    ll
}