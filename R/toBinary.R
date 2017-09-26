#' Transforms a "survival" data frame into a data frame suitable for binary
#' (logistic) regression
#' 
#' The result of the transformation can be used to do survival analysis via
#' logistic regression. If the \code{cloglog} link is used, this corresponds to
#' a discrete time analogue to Cox's proportional hazards model.
#' 
#' toBinary calls \code{risksets} in the \code{eha} package.
#' 
#' @param dat A data frame with three variables representing the survival
#' response. The default is that they are named \code{enter}, \code{exit}, and
#' \code{event}
#' @param surv A character vector with the names of the three variables
#' representing survival.
#' @param strats An eventual stratification variable.
#' @param max.survs Maximal number of survivors per risk set. If set to a
#' (small) number, survivors are sampled from the risk sets.
#' @return Returns a data frame expanded risk set by risk set. The three
#' "survival variables" are replaced by a variable named \code{event} (which
#' overwrites an eventual variable by that name in the input). Two more
#' variables are created, \code{riskset} and \code{orig.row}.
#' \item{event}{Indicates an event in the corresponding risk set.}
#' \item{riskset}{Factor (with levels 1, 2, \ldots{}) indicating risk set.}
#' \item{risktime}{The 'risktime' (age) in the corresponding riskset.}
#' \item{orig.row}{The row number for this item in the original data frame.}
#' @note The survival variables must be three. If you only have \emph{exit} and
#' \emph{event}, create a third containing all zeros.
#' @author Göran Broström
#' @seealso \code{\link[eha]{coxreg}}, \code{\link[stats]{glm}}.
#' @keywords survival cluster
#' @examples
#' 
#' enter <- rep(0, 4)
#' exit <- 1:4
#' event <- rep(1, 4)
#' z <- rep(c(-1, 1), 2)
#' dat <- data.frame(enter, exit, event, z)
#' binDat <- toBinary(dat)
#' dat
#' binDat
#' coxreg(Surv(enter, exit, event) ~ z, method = "ml", data = dat)
#' ## Same as:
#' summary(glm(event ~ z + riskset, data = binDat, family = binomial(link = cloglog)))
#' 
#' @export toBinary
toBinary <- function(dat,
                     surv = c("enter", "exit", "event"),
                     strats,
                     max.survs = NROW(dat)) 
{
    if (!is.data.frame(dat))
      stop("dat must be a data frame")
    if (length(surv) != 3)
      stop("surv must have length 3")
    fixed.names <- names(dat)
    surv.indices <- match(surv, fixed.names)
    if (length(which(is.na(surv.indices)))) {
        x <- which(is.na(surv.indices))
        stop(paste(surv[x], " is not a name in the data frame."))
    }
    enter <- dat[[surv.indices[1]]]
    exit <- dat[[surv.indices[2]]]
    event <- dat[[surv.indices[3]]]

    covars <- dat[, -surv.indices, drop = FALSE]

    nn <- NROW(dat)

    if (missing(strats) || is.null(strats)) strats <- rep(1, nn)
    rs <- risksets(Surv(enter, exit, event), strata = strats, max.survs)

    ## Remove this to keep risksets with no survivors:
    ## (Include this to remove risk sets with no survivors):
    ##weg <- (abs(rs$size - rs$n.events) > 0.01)
    ##rs$riskset <- rs$riskset[rep(weg, rs$size)]
    ##rs$eventset <- rs$eventset[rep(weg, rs$n.events)]
    ##rs$n.events <- rs$n.events[weg]
    ##rs$size <- rs$size[weg]
    ###################
    
    n.rs <- length(rs$size)
    ev <- numeric(sum(rs$size))
    start <- 1
    for (i in 1:n.rs) {
        ev[start:(start + rs$n.events[i] - 1)] <- 1
        start <- start + rs$size[i]
    }
    rs$ev <- ev

    out <- data.frame(event = rs$ev,
                      riskset = factor(rep(1:length(rs$size), rs$size)),
                      ##risktime = rep(rs$risktimes[weg], rs$size)
                      risktime = rep(rs$risktimes, rs$size)
                      )
                      
    out <- cbind(out, covars[rs$riskset, , drop = FALSE])
    out$orig.row <- (1:nn)[rs$riskset]
    out
}
