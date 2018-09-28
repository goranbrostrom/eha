#' Proportional hazards regression with piecewise constant hazards and tabular data.
#' 
#' @param formula a formula
#' @param cuts intervals for constant hazard.
#' @param data a data frame with event, exposure, age plus covariates

#' @export tpchreg
tpchreg <- function(formula, cuts, data){
    ## 'tpch' stands for "Tabular Piecewise Constant Hazard Regression" 
    ## form: event ~ offset(log(exposure)) + age + ....
    ##
    ## cuts: Age intervals
    ##
    ## data: is a data frame containing:
    ##
    ## event: Number of events
    ## exposure: Exposure time
    ## age: Factor of time intervals.
    ## The rest: Covariates
    ## Note:
    ## length(cuts) = length(levels(age)) + 1! (cf: ppch)

    ##form <- as.formula(form)

    if (length(cuts) != length(levels(data$age)) + 1) stop("cuts/age length mismatch")
    fit <- glm(formula, data = data, family = poisson)
    cat(fit$var)
    x <- fit$coefficients
    ages <- which(substr(names(x), 1, 3) == "age")
    n <- length(ages) + 1
    hazards <- numeric(n)
    hazards[1] <- fit$coefficients[1]
    if (n > 1){
        hazards[2:n] <- fit$coefficients[ages] + hazards[1]
    }
    hazards <- exp(hazards)
    ages <- c(1, ages)
    coefficients <- fit$coefficients[-ages]
    out <- list()
    out$coefficients <- coefficients
    out$hazards <- hazards
    out$residuals <- fit$residuals
    out$cuts <- cuts
    out$loglik <- logLik(fit)
    out$terms <- fit$terms
    out$aic <- extractAIC(fit)
    out$call <- fit$call
    out$pfixed = FALSE
    class(out) <- "tpchreg"
    out$call <- match.call()
    out
}

#' @export
extractAIC.tpchreg <- function(x, scale, k = 2, ...){
    x$aic
}

#' @export
nobs.tpchreg <- function(x, ...){
    length(x$residuals)
}
