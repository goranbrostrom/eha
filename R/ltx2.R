#' LaTeX alternative printing of regression results.
#' 
#' This (generic) function prints the LaTeX code of the results of a fit from
#' \code{\link{coxreg}}, \code{\link{phreg}}, \code{\link{tpchreg}}, 
#' or \code{\link{aftreg}}.
#' 
#' @param x The output from a call to \code{coxreg}, \code{tpchreg}, or
#' \code{aftreg}
#' @param caption A suitable caption for the table.
#' @param label A label used in the LaTeX code.
#' @param dr Output from a \code{drop1} call.
#' @param digits Number of digits to be printed.
#' @param conf Confidence intervals level.
#' @param keep Number of covariates to present.
#' @param \dots Not used.
#' @return LaTeX code version of the results from a run with
#' \code{\link{coxreg}}, \code{\link{phreg}}, \code{\link{phreg}}, 
#' \code{\link{aftreg}}.
#' @note Resulting tables contain estimated hazard ratios and confidence limits
#' instead of regression coefficients and standard errors as in \code{\link{ltx}}.
#' @author Göran Broström.
#' @seealso \code{xtable}, \code{\link{coxreg}}, \code{\link{phreg}}, 
#' \code{\link{phreg}}, \code{\link{aftreg}}, and \code{\link{ltx}}.
#' @keywords printing
#' @examples
#' 
#' data(oldmort)
#' fit <- coxreg(Surv(enter, exit, event) ~ sex, data = oldmort)
#' ltx2(fit, caption = "A test example.", label = "tab:test1") 
#' 
#' @export ltx2
ltx2 <- function(x,
                 caption = NULL,
                 label = NULL,
                 dr = NULL,
                 digits = max(options()$digits - 4, 4),
                 conf = 0.95,
                 keep = NULL, ...){
    UseMethod("ltx2")
}

#' @export
ltx2.coxreg <- function(x, caption = NULL, label = NULL, dr = NULL, 
          digits=max(options()$digits - 4, 3), conf = 0.95, keep = NULL, ...){

    ## keep is checked in 
    if (is.null(dr)){
         dr <- drop1(x, test = "Chisq")
    }
    if (!inherits(x, "coxreg")){
        stop("only for coxreg objects")
    }
    if (!is.null(cl<- x$call)) {
	##cat("Call:\n")
	##dput(cl)
	##cat("\n")
	}

    if (!is.null(x$fail)) {
        if (x$fail != 0){
            cat(" coxreg failed with: ")
            stop(x$fail)
        }
    }

    if (!length(x$coefficients)){
        cat("Null log likelihood = $", x$loglik[2], "$\n")
        return()
    }

    savedig <- options(digits = digits)
    on.exit(options(savedig))

    ltxCoef3(x, dr, conf, keep, digits, caption) # Print regression coefficients
### New 2020-11-26:    
    if (inherits(x, "summary.tpchreg")){
        ivl <- paste("(", min(x$cuts), ", ", max(x$cuts), "]", sep = "")
        if (x$n.strata == 1){
            cat("Restr. mean survival: & ", x$rmean, " & in", ivl, "\\\\ \n")
        }else{
            names(x$rmean) <- x$strata
            cat("Restr. mean survival in & ", ivl, ": \\\\ \n")
            for (ll in 1:(x$n.strata - 1)){
                cat(" & ", formatC(x$strata[ll], format = "g", flag = "#"))
            }
            cat(" & ", x$strata[x$n.strata], "\\\\ \n")
            for (ll in 1:(x$n.strata - 1)){
                cat(" & ", formatC(x$rmean[ll], format = "g", flag = "#"))
            }
            cat(" & ", formatC(x$rmean[x$n.strata], format = "g", flag = "#"), "\\\\ \n")
        }
    }
### End New 2020-11-26. 
    if (TRUE){
    cat("\\hline \n")

    cat("\\end{tabular}\n")

    if (!is.null(label)){
        cat("\\label{", label, "} \n", sep = "")
    }
    cat("\\normalsize \n")
    cat("\\end{center} \n")
    cat("\\end{table} \n\n\n")
    cat(" \n\n")
    } # End FALSE
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
    if (FALSE){
        logtest <- -2 * (x$loglik[1] - x$loglik[2])
        if (is.null(x$df)) df <- sum(!is.na(coef))
        else  df <- round(sum(x$df),2)
        cat("\n")
        cat(formatC("Events", width = 25, flag = "-"), x$events, "\n")
        cat(formatC("Total time at risk", width = 25, flag = "-"),
            formatC(x$ttr, digits = 5, format = "fg"), "\n")
        cat(formatC("Max. log. likelihood", width = 25, flag = "-"),
            formatC(x$loglik[2], digits = 5, format = "fg"), "\n")
        cat(formatC("LR test statistic", width = 25, flag = "-"),
            format(round(logtest, 2)), "\n")
        cat(formatC("Degrees of freedom", width = 25, flag = "-"),
            formatC(df, digits = 0, format = "f"), "\n")
        cat(formatC("Overall p-value", width = 25, flag = "-"),
            format.pval(1 - pchisq(logtest, df), digits = 6, "\n"))
        cat("\n")
        if (length(x$icc))
            cat("   number of clusters=", x$icc[1],
                "    ICC=", format(x$icc[2:3]), "\n")
        invisible(x)
    }
}

#' @export
ltx2.phreg <- function(x, caption = NULL, label = NULL, dr = NULL, 
          digits=max(options()$digits - 4, 3), conf, keep = NULL, ...){
    
    if (!("phreg" %in% class(x))) stop("Only for 'phreg' output")
    if (is.null(dr)){
        dr <- drop1(x, test = "Chisq")
    }

    if (!is.null(cl<- x$call)) {
	##cat("Call:\n")
	##dput(cl)
	##cat("\n")
	}

    if (!is.null(x$fail)) {
        if (x$fail != 0){
            cat(" phreg failed with: ")
            stop(x$fail)
        }
    }

    if (!length(x$coefficients)){
        cat("Null log likelihood = ", x$loglik[2], "\n")
        return()
    }

    if (x$pfixed){

        n.slsh <- 1

    }else{
        n.slsh <- 2 * x$n.strata

    }

    savedig <- options(digits = digits)
    on.exit(options(savedig))

    ltxCoef3(x, dr, conf, keep, digits, caption) # Print regression coefficients.


    cat("\\end{tabular}\n")

    if (!is.null(label)){
        cat("\\label{", label, "} \n", sep = "")
    }

    cat("\\end{center} \n")
    cat("\\end{table} \n\n\n")
    cat(" \n\n")
#####################################

    if (!is.null(x$frailty)){
        cat("\nFrailty standard deviation = ", x$sigma, "\n")
        cat("                      S.E. = ", x$sigma.sd, "\n\n")
    }

}

ltxCoef3 <- function(x, dr, conf, keep, digits, caption){
    coef <- x$coef

    se <- sqrt(diag(x$var))

##    wald.p <- formatC(1 - pchisq((coef/ se)^2, 1),
##                      digits = digits,
##                      format = "f")
    if(is.null(coef) | is.null(se))
        stop("Input is not valid")
    ## Check for dr:
    lp <- TRUE
    if (is.null(dr)) {
        dr <- drop1(x, test = "Chisq")
    }
    if (lp){
        lpval <- formatC(dr[-1, 4], digits = digits, format = "f")
    }
    
    qn <- qnorm(1 - (1 - conf) / 2)
    lower <- exp(coef - qn * se)
    upper <- exp(coef + qn * se)
    
#####################################
    cat("\\begin{table}[ht] \n")
    if (!is.null(caption)){
        cat("\\caption{", caption, "} \n", sep = "")
    }
    
    
    cat("\\begin{center} \n")
    cat("\\footnotesize \n") # NOTE!!
    cat("\\begin{tabular}{lrrrrr} \n")
    cat("\\hline \n")
    if (lp){
        if ("aftreg" %in% x$class){
            if (x$param == "default"){
                cat("\\bf Covariate & \\bf Mean & \\bf Coef & \\bf Time accn. & \\bf S.E. & \\bf  L-R p \\\\ \\hline\n")
            }else{
                cat("\\bf Covariate & \\bf Mean & \\bf Coef & \\bf Life expn. & \\bf S.E. &   \\bf L-R p \\\\ \\hline\n")
            }
        }else{
            cat("\\bf Covariate & \\bf Mean & \\bf H.R. & \\bf lowCI &  \\bf highCI & \\bf L-R p \\\\ \\hline\n")
        }
    }else{
        if ("aftreg" %in% x$class){
            if (x$param == "default"){
                cat("\\bf Covariate & \\bf Mean & \\bf Coef & \\bf Time accn. & \\bf S.E. &  \\bf Wald p \\\\ \\hline\n")
            }else{
                cat("\\bf Covariate & \\bf Mean & \\bf Coef & \\bf Life expn. & \\bf S.E. &  \\bf Wald p \\\\ \\hline\n")
            }
        }else{
            cat("\\bf Covariate & \\bf Mean & \\bf Coef & \\bf H.R. & \\bf S.E. &   \\bf Wald p \\\\ \\hline\n")
        }
    }
    e.coef <- formatC(exp(coef), digits = digits, format = "f")
    lower <- formatC(lower, digits = digits, format = "f")
    upper <- formatC(upper, digits = digits, format = "f")
    
    ett <-  1L
##    noll <-  0L

    factors <- attr(x$terms, "factors")
    resp <- attr(x$terms, "response")
    row.strata <- attr(x$terms, "specials")$strata
    if (!is.null(row.strata))
      col.strata <- which(factors[row.strata, ] == 1)
    else col.strata <- NULL
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

    ord <- attr(x$terms, "order")
    if (!is.null(col.strata)) ord <- ord[-col.strata]

    index <- 0
    lpindx <- 0

    if (is.null(keep) || missing(keep)){
        goto <- length(term.names)
    }else{
        goto <- min(keep, length(term.names))
    }
    
    for (term.no in 1:goto){

        if (ord[term.no] == 1){
            covar.no <- which(factors[, term.no] == 1)

            if (isF[covar.no]){ ## Factors:
                if (lp){
                    lpindx <- lpindx + 1
                    cat(covar.names[covar.no], " & & & & &",
                        lpval[lpindx], " \\\\ \n")
                }else{
                    cat(covar.names[covar.no], " \\\\ \n")
                }
                ##p <- match(covar.names[covar.no], names(data))
                no.lev <- length(x$levels[[covar.no]])
                x$levels[[covar.no]] <-
                    substring(x$levels[[covar.no]], 1, 16)
                lb <- paste("\\em", x$levels[[covar.no]][1], sep = " ")
                cat("\\multicolumn{1}{r}{", lb, "}",
                ## cat(formatC(x$levels[[covar.no]][1], width = 16, flag = "+"),
                    " & ",
                    formatC(x$w.means[[covar.no]][1],
                            width = 8, digits = digits, format = "f"), " & ",
                    
                    ett, " & ",
                    "\\multicolumn{2}{c}{(reference)} \\\\ \n", sep = "")
                for (lev in 2:no.lev){
            ##cat("lev = ", lev, "\n")
                    index <- index + 1
                    lb <- paste("\\em", x$levels[[covar.no]][lev], sep = " ")
                    cat("\\multicolumn{1}{r}{", lb, "}",
                    ##cat(formatC(x$levels[[covar.no]][lev], width = 16,
                      ##          flag = "+"),
                        " & $",
                        formatC(x$w.means[[covar.no]][lev],
                                width = 8, digits = digits, format = "f"), "$ & $",
                        ##coef[index], "$ & $",
                        e.coef[index], "$ & $",
                        lower[index], "$ & $", 
                        upper[index], "$") 
                        ##formatC(" ", width = 9),
                    
                        cat("\\\\ \n")
                   
                }
            }else{ ## Covariates:
                index <- index + 1
                cat(covar.names[covar.no], 
                ##cat(formatC(substr(covar.names[covar.no], 1, 16),
                  ##          width = 16, flag = "-"),
                    " & $ ",
                    formatC(x$w.means[[covar.no]],
                            width = 8, digits = digits, format = "f"),
                    " $ & $ ",
                    ##coef[index], " $ & $ ",
                    e.coef[index], " $ & $ ",
                                        #exp(coef[index]),
                    lower[index], "$ & $",
                    upper[index], " $ & $ ")
                if (lp){
                    lpindx <- lpindx + 1
                    ppv <- lpval[lpindx]
                }
                cat(ppv, "$ \\\\ \n")
            }
        }else if (ord[term.no] > 1){ ## Interactions:
            cat(formatC(term.names[term.no], width = 16, flag = "-"), " \\\\ \n")
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
                cat(formatC(" ", width = 2),
                    formatC(substring(vn, 1, 22), width = 22, flag = "-"),
                    ##format(" ", 8),
                    " &  & $ ",
                    ##coef[index], " & ", 
                    e.coef[index], "$ & $",
                    lower[index], "$ & $",
                    upper[index], "$ & $",
                    formatC(" ", width = 9),
                    ##formatC(wald.p[index],
                
                      ##      digits = digits,
                        ##    width = digits + 2,
                          ##  format = "f"),
                    ##signif(1 - pchisq((coef[index]/ se[index])^2, 1), digits - 1),
                    "$ \\\\ \n")
            }
        }
        
    }
    cat("\\hline \n")
    cat("Events & ", x$n.events, " & TTR & ", x$ttr, "\\\\ \n")
    cat("Max. logLik. & $ ", x$loglik[2], " $ & ", "Conf level & ",   conf,  "\\\\ \\hline \n")
    
    cat("\\hline \n")
    ## Remove the rest?
##    cat("\\end{tabular}\n")
##    if (!is.null(label)) {
##        cat("\\label{", label, "} \n", sep = "")
##    }
##    cat("\\normalsize \n")
##    cat("\\end{center} \n")
##    cat("\\end{table} \n\n\n")
##    cat(" \n\n")
}
