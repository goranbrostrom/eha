#' Retrieves regression tables
#' 
#' @param x A \code{summary.XXreg} object, typically the result of running
#' \code{summary.XXreg}, summary on a XXreg object.
#' @param digits Output format.
#' @param short If TRUE, return only coefficients table.
#' @param \dots Other arguments.
#' @return A character data frame, ready to print in various formats.
#' @note Should not be used if interactions present.
#' @author Göran Broström
#' @seealso \code{\link{coxreg}}, \code{\link{summary.coxreg}}
#' @keywords survival
#' @export
regtable <- function(x, digits = 3, short = TRUE, ...){

    if (!("summary.coxreg" %in% class(x))){
    ##    stop("Only for 'summary.coxreg' ojects.")
    }

    ##if (!is.null(cl<- x$call)) {
	##cat("Call:\n")
	##dput(cl)
	##cat("\n")
	##}

    if (!is.null(x$fail)) {
        if (x$fail != 0){
            cat(" coxreg failed with: ")
            stop(x$fail)
        }
    }

    if (x$nullModel){
        ##cat("Null log likelihood = ", x$loglik[2], "\n")
        cat("Null model\n")
        return()
    }

    savedig <- options(digits = digits)
    on.exit(options(savedig))

    coef <- x$coefficients[, 1]
    se <- x$coefficients[, 3]

    wald.p <- formatC(pchisq((coef / se)^2, 1, lower.tail = FALSE),
                      digits = digits,
                      width = 9, format = "f")
    if(is.null(coef) | is.null(se))
        stop("Input is not valid")
    ## Check for dr:
    lp <- !is.null(x$dr)
    if (lp){
        dr <- x$dr[rownames(x$dr) %in% x$covars, ]
        lpval <- formatC(dr[, 4], digits = digits, width = 9, format = "f")
    }
    ord <- attr(x$terms, "order")
    ## New, 2020-07-26:
    if (any(ord > 1)){
        lp <- FALSE
        if (!is.null(x$dr)){
            print(x$dr)
            cat("\n")
        }
        return(NULL)
    }
    ####
    
    res <- matrix(NA, ncol = 7, nrow = 0)
    colnames(res) <- c("Covariate", "level", "W_mean", "Coef", "HR", "SE", "LR_p")
    if (inherits(x, "summary.aftreg")){
        if (x$param == "lifeAcc"){
            colnames(res)[5] <- "lifeAcc"
        }else{
            colnames(res)[5] <- "lifeExp"
        }
    }
    if (!lp) colnames(res)[7] <- "Wald_p"
    ##res <- as.data.frame(res)
    ##return(res)
#####################################
    e.coef <- formatC(exp(coef), width = 9, digits = digits, format = "f")
    coef <- formatC(coef, width = 9, digits = digits, format = "f")
    se <- formatC(se, width = 9, digits = digits, format = "f")
    
    ett <-  formatC(1, width = 9, digits = 0, format = "f")
    noll <-  formatC(0, width = 5, digits = 0, format = "f")

    factors <- attr(x$terms, "factors")
    resp <- attr(x$terms, "response")
    row.strata <- attr(x$terms, "specials")$strata
    if (!is.null(row.strata)){
         col.strata <- which(factors[row.strata, ] == 1)
    }else{
        col.strata <- NULL
    }
    if (!is.null(col.strata)){
        factors <-
            attr(x$terms, "factors")[-c(resp, row.strata), -col.strata,
                                 drop = FALSE]
    }else{
        factors <-
            attr(x$terms, "factors")[-c(resp, row.strata), ,
                                     drop = FALSE]
    }

    covar.names <- x$covars
    term.names <- colnames(factors)

    isF <- x$isF

    
    if (!is.null(col.strata)) ord <- ord[-col.strata]

    index <- 0
    lpindx <- 0
    row <- 0
    for (term.no in 1:length(term.names)){

        if (ord[term.no] == 1){
            covar.no <- which(factors[, term.no] == 1)

            if (isF[covar.no]){ ## Factors:
                if (lp){
                    lpindx <- lpindx + 1
                    outp <- formatC(covar.names[covar.no], width = 56, 
                                    flag = "-")
                    row <- row + 1
                    res <- rbind(res, rep("", 7))
                    res[row, 1] <- covar.names[covar.no]
                    res[row, 7] <- lpval[lpindx]
                    ##return(res)
                    ##cat(outp, lpval[lpindx], "\n")
                }else{
                    row <- row + 1
                    res <- rbind(res, rep("", 7))
                    res[row, 1] <- covar.names[covar.no]
                    ##cat(covar.names[covar.no], "\n")
                }
                ##p <- match(covar.names[covar.no], names(data))
                no.lev <- length(x$levels[[covar.no]])
                x$levels[[covar.no]] <-
                    substring(x$levels[[covar.no]], 1, 16)
                row <- row + 1
                res <- rbind(res, rep("", 7))
                res[row, 2] <- x$levels[[covar.no]][1]
                res[row, 3] <- formatC(x$w.means[[covar.no]][1],
                                       width = 10, digits = 3, format = "f")
                res[row, 4] <- noll
                res[row, 5] <- ett
                ##cat(formatC(x$levels[[covar.no]][1], width = 16, flag = "+"),
                  ##  formatC(x$w.means[[covar.no]][1],
                    ##        width = 10, digits = 3, format = "f"),
                ##    noll,
                 ##   ett,
                   ## "(reference)\n")
                for (lev in 2:no.lev){
            ##cat("lev = ", lev, "\n")
                    index <- index + 1
                    row <- row + 1
                    res <- rbind(res, rep("", 7))
                    res[row, 2] <- x$levels[[covar.no]][lev]
                    res[row, 3]  <- formatC(x$w.means[[covar.no]][lev], width = 10,
                                            digits = 3, format = "f")
                    res[row, 4] <- coef[index]
                    res[row, 5] <- e.coef[index]
                    res[row, 6] <- se[index]
                    if (!lp) res[row, 7] <- wald.p[index]
                    if (FALSE){
                    cat(formatC(x$levels[[covar.no]][lev], width = 16,
                                flag = "+"),
                        formatC(x$w.means[[covar.no]][lev],
                                width = 10, digits = 3, format = "f"),
                        coef[index],
                        e.coef[index],
                        se[index])
                    if (lp){
                        cat("\n")
                    }else{
                        ##formatC(" ", width = 9),
                        cat(formatC(wald.p[index],
                                digits = 3,
                                width = digits + 2,
                                format = "f"),
                        ##signif(1 - pchisq((coef/ se)^2, 1), digits - 1),
                        "\n")
                    }
                    }
                }
            }else{ ## Covariates:
                index <- index + 1
                xxx <- x$w.means[[covar.no]]
                if (inherits(xxx, "Date")){
                    xxx <- as.character(xxx)
                }else{
                    xxx <- formatC(xxx,
                                   width = 10, digits = 3, format = "f")
                }
                row <- row + 1
                res <- rbind(res, rep("", 7))
                res[row, 1] <- covar.names[covar.no]
                res[row, 3] <- xxx
                res[row, 4] <- coef[index]
                res[row, 5] <- e.coef[index]
                res[row, 6] <- se[index]
                ##cat(formatC(substr(covar.names[covar.no], 1, 16),
                  ##          width = 16, flag = "-"),
                    ##xxx,
                ##    coef[index],
                ##    e.coef[index],
                                        #exp(coef[index]),
                  ##  se[index])
                ##formatC(" ", width = 9),
                if (lp){
                    lpindx <- lpindx + 1
                    ppv <- lpval[lpindx]
                }else{
                    ppv <- wald.p[index]
                }
                res[row, 7] <- ppv
                ##cat(formatC(ppv,
                  ##          digits = 3,
                    ##        width = digits + 2,
                      ##      format = "f"), "\n")
            }
        }else if (ord[term.no] > 1){ ## Interactions:
            row <- row + 1
            rbind(res, rep("", 7))
            res[row, 1] <- term.names[term.no]
            ##cat(formatC(term.names[term.no], width = 16, flag = "-"), "\n")
            niv <- numeric(ord[term.no])
            covar.no <- which(factors[, term.no] == 1)

            for (i in 1:ord[term.no]){
                if (isF[covar.no[i]]){
                    niv[i] <- length(x$levels[[covar.no[i]]]) - 1
                }else{
                    niv[i] <- 1
                }
            }
            stt <- index + 1
            for (index in stt:(stt + prod(niv) - 1)){
                vn <- sub(covar.names[covar.no[1]], "", names(coef)[index])
                for (i in 1:ord[term.no]){
                    vn <- sub(covar.names[covar.no[i]], "", vn)
                }
                ##          cat(format(names(coef)[index], 15, "+"),
                row <- row + 1
                res <- rbind(res, rep("", 7))
                res[row, 1] <- substring(vn, 1, 22)
                res[row, 4] <- coef[index]
                res[row, 5] <- e.coef[index]
                res[row, 6] <- se[index]
                if (!lp) res[row, 7] <- wald.p[index]
                if (FALSE){
                cat(formatC(" ", width = 2),
                    formatC(substring(vn, 1, 22), width = 22, flag = "-"),
                    ##format(" ", 8),
                    coef[index],
                    e.coef[index],
                    se[index],
                    ##formatC(" ", width = 9),
                    formatC(wald.p[index],
                            digits = 3,
                            width = digits + 2,
                            format = "f"),
                    ##signif(1 - pchisq((coef[index]/ se[index])^2, 1), digits - 1),
                    "\n")
                }
            }
        }
        
    }
      
#####################################
    if(FALSE){
        tmp <- cbind(coef, exp(coef), se,
                     signif(1 - pchisq((coef/ se)^2, 1), digits - 1))
        dimnames(tmp) <- list(names(coef), c("coef", "rel. risk",
                                             "se(coef)", "p"))
        
        cat("\n")
        prmatrix(tmp)
    }

    if (!is.null(x$frailty)){
        cat("\nFrailty standard deviation = ", x$sigma, "\n")
        cat("                      S.E. = ", x$sigma.sd, "\n\n")
    }
    
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    if (is.null(x$df)){
        df <- sum(!is.na(coef))
    }else{
        df <- round(sum(x$df),2)
    }
    if (FALSE){ #if (!short){
        cat("\n")
        cat(formatC("Events", width = 25, flag = "-"), x$n.events, "\n")
        cat(formatC("Total time at risk", width = 25, flag = "-"),
            formatC(x$ttr, digits = 5, format = "fg"), "\n")
        cat(formatC("Max. log. likelihood", width = 25, flag = "-"),
            formatC(x$loglik[2], digits = 5, format = "fg"), "\n")
        cat(formatC("LR test statistic", width = 25, flag = "-"),
            format(round(logtest, 2), nsmall = 2), "\n")
        cat(formatC("Degrees of freedom", width = 25, flag = "-"),
            formatC(df, digits = 0, format = "f"), "\n")
        cat(formatC("Overall p-value", width = 25, flag = "-"),
            format.pval(1 - pchisq(logtest, df), digits = 6, "\n"))
        cat("\n")
        if (length(x$icc))
            cat("   number of clusters=", x$icc[1],
                "    ICC=", format(x$icc[2:3]), "\n")
    }
    as.data.frame(res) ##invisible(x)
}