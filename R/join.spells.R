#' Straighten up a survival data frame
#' 
#' Unnecessary cut spells are glued together, overlapping spells are
#' "polished", etc.
#' 
#' In case of overlapping intervals (i.e., a data error), the appropriate id's
#' are returned if \code{strict} is \code{TRUE}.
#' 
#' @param dat A data frame with names enter, exit, event, id.
#' @param strict If TRUE, nothing is changed if errors in spells (non-positive
#' length, overlapping intervals, etc.) are detected. Otherwise (the default),
#' bad spells are removed, with "earlier life" having higher priority.
#' @param eps Tolerance for equality of two event times. Should be kept small.
#' @return A data frame with the same variables as the input, but individual
#' spells are joined, if possible (identical covariate values, and adjacent
#' time intervals).
#' @author Göran Broström
#' @seealso \code{\link{coxreg}}, \code{\link{aftreg}},
#' \code{\link{check.surv}}
#' @references Therneau, T.M. and Grambsch, P.M. (2000). \emph{Modeling
#' Survival Data: Extending the Cox model.} Springer.
#' @keywords manip survival
#' @export join.spells
join.spells <- function(dat, strict = FALSE, eps = 1.e-8){
    ## Survival data: (enter, exit], event (0-1, or TRUE/FALSE),
    ## birthdate in years since 1 jan 0, eg 1877.500 = 1 july 1877
    ## Assumes: enter, exit, event, id, birthdate
    ## Must have unique id (as a covariate).

    
    resp <- match(c("enter", "exit", "event"), names(dat))
    if (any(is.na(resp))) stop("Wrong variable names")
    
    koll <- match(c("id"), names(dat))
    if (any(is.na(resp))) stop("No 'id' in variable names")

    ## First, if strict, check data:
    if (strict){
        res.check <- check.surv(dat$enter, dat$exit, dat$event, dat$id)
        if (length(res.check)){
            cat("Error in individual(s). Return value is id of the bad.\n")
            return(res.check)
        }
    }
    covar <- dat[ , -resp]
    if (na.cov <- any(is.na(covar))){
         warning("NA's in covariates: temporarily replaced by -999")
        covar[is.na(covar)] <- -999
    }
    n.cov <- ncol(covar)
    n.rec <- nrow(covar)
    all <- unique(dat$id)
    nn <- length(all)
    
    ide <- as.integer(factor(dat$id, labels = 1:nn))
    ord <- order(ide, dat$enter, dat$exit)
    dat <- dat[ord, ]
    
    res <- .Fortran("cleanup",
                    as.double(t(covar)),
                    as.double(dat$enter),
                    as.double(dat$exit),
                    as.integer(dat$event),
                    as.integer(ide),
                    as.integer(n.cov),
                    as.integer(n.rec),
                    as.integer(nn),
                    as.double(eps),
                    new.n.rec = integer(1),
                    new.cov = double(n.rec * n.cov),
                    enter = double(n.rec),
                    exit = double(n.rec),
                    event = integer(n.rec),
                    id = integer(n.rec)
                    ##DUP = FALSE,
                    )
    
  
    out <- data.frame(new.id = res$id[1:res$new.n.rec],
                      enter = res$enter[1:res$new.n.rec],
                      exit = res$exit[1:res$new.n.rec],
                      event = res$event[1:res$new.n.rec]
                      )
    
    new.cov <-
        data.frame(matrix(res$new.cov, byrow = TRUE,
                          ncol = n.cov))[1:res$new.n.rec, , drop = FALSE]
    
    names(new.cov) <- names(covar)
    if (na.cov) covar[covar == -999] <- NA
    cbind(out, new.cov)
}
