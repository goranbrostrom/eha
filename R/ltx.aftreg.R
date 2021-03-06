#' @export
ltx.aftreg <- function(x, caption = NULL, label = NULL, dr = NULL, 
          digits=max(options()$digits - 4, 3), ...){
    
    if (!("aftreg" %in% class(x))) stop("Only for 'aftreg' output")

    if (!is.null(cl<- x$call)) {
	##cat("Call:\n")
	##dput(cl)
	##cat("\n")
	}

    if (!is.null(x$fail)) {
        if (x$fail != 0){
            cat(" aftreg failed with: ")
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

    coef <- x$coefficients
    se <- sqrt(diag(x$var))
    
    wald.p <- formatC(1 - pchisq((coef/ se)^2, 1),
                      digits = digits,
                      format = "f")
    if(is.null(coef) | is.null(se))
        stop("Input is not valid")
    ## Check for dr:
    lp <- !is.null(dr) 
    if (lp){
        lpval <- formatC(dr[-1, 4], digits = digits, format = "f")
    }
#####################################
    cat("\\begin{table}[ht] \n")
    cat("\\begin{center} \n")
    cat("\\footnotesize \n")
    if (!is.null(caption)){
        cat("\\caption{", caption, "} \n", sep = "")
    }
    
    if (!is.null(label)){
        cat("\\label{", label, "} \n", sep = "")
    }
    
    cat("\\begin{tabular}{lrrrrr} \n")
    cat("\\hline \n")
    if (lp){
        ptxt <- "\\bf L-R $p$"
    }else{
        ptxt <- "\\bf Wald $p$"
    }
    if(x$param == "lifeExp"){
        etxt <- "\\bf Ext'd life"
    }else{
        etxt <- "\\bf Acc'd time"
    }
    if (!is.null(x$covars)){
        cat(paste("\\bf Covariate & \\bf Mean & \\bf Coef &", 
                  etxt, " & \\bf S.E. & ", ptxt,
                  " \\\\ \\hline\n"))
    }

    e.coef <- formatC(exp(coef), digits = digits, format = "f")
    coef <- formatC(coef, digits = digits, format = "f")
    se <- formatC(se, digits = digits, format = "f")
    
    ett <-  1L
    noll <-  0L
    
    factors <- attr(x$terms, "factors")
    resp <- attr(x$terms, "response")
    row.strata <- attr(x$terms, "specials")$strata
    if (!is.null(row.strata))
      col.strata <- which(factors[row.strata, ] == 1)
    else col.strata <- NULL
    if (!is.null(x$covars)){
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
     }   
    index <- 0
    lpindx <- 0
    
    if (FALSE){ ## Do NOT introduce (Intercept) for now!
        cat("(Intercept)", 
            ##cat(formatC(substr(covar.names[covar.no], 1, 16),
            ##          width = 16, flag = "-"),
            " & ",
            " ",
            " & ",
            coef[index], " & ",
            e.coef[index], " & ",
                                        #exp(coef[index]),
            se[index], " & ")
        if (lp){
            lpindx <- lpindx + 1
            ppv <- lpval[lpindx]
        }else{
            ppv <- wald.p[index]
        }
        ##formatC(" ", width = 9),
        cat(ppv,
            ##signif(1 - pchisq((coef/ se)^2, 1), digits - 1),
            "\\\\ \n")
    }
    if (!is.null(x$covars)){
        for (term.no in 1:length(term.names)){
            
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
                        " & $ ",
                        formatC(x$w.means[[covar.no]][1],
                                width = 8, digits = digits, format = "f"), " $ & ",
                        noll, " & ",
                        ett, " & ",
                        "\\multicolumn{2}{c}{(reference)} \\\\ \n", sep = "")
                    for (lev in 2:no.lev){
                        ##cat("lev = ", lev, "\n")
                        index <- index + 1
                        lb <- paste("\\em", x$levels[[covar.no]][lev], sep = " ")
                        cat("\\multicolumn{1}{r}{", lb, "}",
                            ##cat(formatC(x$levels[[covar.no]][lev], width = 16,
                            ##          flag = "+"),
                            " & $ ",
                            formatC(x$w.means[[covar.no]][lev],
                                    width = 8, digits = digits, format = "f"), " $ & $ ",
                            coef[index], " $ & $ ",
                            e.coef[index], " $ & $ ",
                            se[index]) 
                        ##formatC(" ", width = 9),
                        if (lp){
                            cat(" $ \\\\ \n")
                        }else{
                            cat(" $ & $ ", wald.p[index],
                                " $ \\\\ \n", sep = "")
                        }
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
                        coef[index], " $ & $ ",
                        e.coef[index], " $ & $ ",
                                        #exp(coef[index]),
                        se[index], " $ & $ ")
                    if (lp){
                        lpindx <- lpindx + 1
                        ppv <- lpval[lpindx]
                    }else{
                        ppv <- wald.p[index]
                    }
                    ##formatC(" ", width = 9),
                    cat(ppv,
                        ##signif(1 - pchisq((coef/ se)^2, 1), digits - 1),
                        " $ \\\\ \n")
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
                        " &  & ",
                        coef[index], " & ", 
                        e.coef[index], " & ",
                        se[index], " & ",
                        ##formatC(" ", width = 9),
                        formatC(wald.p[index],
                                digits = digits,
                            width = digits + 2,
                                format = "f"),
                        ##signif(1 - pchisq((coef[index]/ se[index])^2, 1), digits - 1),
                        "\\\\ \n")
                }
            }
            
        }
    }## if (!is.null(x$covars)
    cat("\\hline \n")
    cat("Baseline parameters \\\\\n")
    if (FALSE){ ## Skip scale/shape printing ...
        for (i in 1:n.slsh){
            jup <- length(coef)
            ss.names <- names(coef[(jup - n.slsh + 1):jup])
            index <- index + 1
            ## covar.no <- covar.no + 1
            cat(formatC(ss.names[i], width = 16, flag = "-"),
                formatC(" & & ",
                        width = 8, digits = 3, format = "c"),
                coef[index], " & ",
                e.coef[index], " & ",
                #exp(coef[index]),
                se[index], " & ",
                #formatC(" ", width = 1),
                formatC(wald.p[index],
                        digits = 3,
                        width = digits + 2,
                        format = "f"), " \\\\ ", 
                ##signif(1 - pchisq((coef/ se)^2, 1), digits - 1),
                "\n")
        }
    }
    if (x$n.strata == 1){
        cat("Baseline expectation: & ", x$baselineMean, "\\\\", "\\hline", "\n")
    }else{
        for (j in 1:x$n.strata){
            cat("Baseline expectation", j, ": &", x$baselineMean[j], "\\\\", "\n")
        }
        cat("\\hline \n")
    }
    cat("\\hline \n")
    cat("Events & ", x$n.events, " & TTR & ", x$ttr, "\\\\ \n")
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    df <- x$df
    if (is.null(df)){
        pvale <- 1
    }else{
        pvale <- pchisq(logtest, df, lower.tail = FALSE)
        pvale <- round(pvale, digits = 6)
    }
    cat("Max. Log Likelihood & $ ", x$loglik[2], " $ & $p$-value & ", pvale,
        "\\\\ \\hline \n")
    cat("\\hline \n")

    cat("\\end{tabular}\n")


    cat("\\end{center} \n")
    cat("\\end{table} \n\n\n")
    cat(" \n\n")
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
