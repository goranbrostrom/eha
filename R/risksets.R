#' Finds the compositions and sizes of risk sets
#' 
#' Focus is on the risk set composition just prior to a failure.
#' 
#' If the input argument max.survs is left alone, all survivors are accounted
#' for in all risk sets.
#' 
#' @param x A \code{Surv} object.
#' @param strata Stratum indicator.
#' @param max.survs Maximum number of survivors in each risk set. If smaller
#' than the 'natural number', survivors are sampled from the present ones. No
#' sampling if missing.
#' @param members If TRUE, all members of all risk sets are listed in the
#' resulting list, see below.
#' @param collate_sets logical. If TRUE, group information by
#' risk sets in a list. Only if \code{members = TRUE}.
#' @return A list with components (if \code{collate_sets = FALSE})
#' \item{antrs}{No. of risk sets in each
#' stratum. The number of strata is given by \code{length(antrs)}.}
#' \item{risktimes}{Ordered distinct failure time points.}
#' \item{eventset}{If
#' 'members' is TRUE, a vector of pointers to events in each risk set, else
#' NULL.}
#' \item{riskset}{If 'members' is TRUE, a vector of pointers to the
#' members of the risk sets, in order. The 'n.events' first are the events. If
#' 'members' is FALSE, 'riskset' is NULL.}
#' \item{size}{The sizes of the risk
#' sets.}
#' \item{n.events}{The number of events in each risk set.}
#' \item{sample_fraction}{If 'members' is TRUE, the sampling fraction of
#' survivors in each risk set.
#' }
#' @note Can be used to "sample the risk sets". 
#' @author Göran Broström
#' @seealso \code{\link{table.events}}, \code{\link{coxreg}}.
#' @keywords survival
#' @examples
#' 
#'  enter <- c(0, 1, 0, 0)
#'  exit <- c(1, 2, 3, 4)
#'  event <- c(1, 1, 1, 0)
#'  risksets(Surv(enter, exit, event))
#' 
#' @export risksets
#' 

risksets <- function (x, strata = NULL, max.survs = NULL, members = TRUE,
                      collate_sets = FALSE){
    ## x is a Surv (survival) object.

  nn <- NROW(x)
  if (is.null(strata)){
      strata <- rep(1, nn)
  }else{
      if (length(strata) != nn) stop("'strata' has wrong length")
      else
          strata <- as.integer(factor(strata))
  }

  if (is.null(max.survs)) max.survs <- nn - 1
  if (NCOL(x) == 2){
      enter <- numeric(nn)
      exit <- x[, 1]
      event <- (x[, 2] != 0)
  }else{
      if (NCOL(x) != 3) stop("'x' is not a Surv object")
      enter <- x[, 1]
      exit <- x[, 2]
      event <- (x[, 3] != 0) ## TRUE == event
  }

  ord <- order(strata, exit, -event)
  strata <- strata[ord]
  enter <- enter[ord]
  exit <- exit[ord]
  event <- event[ord]
  w.totrs <- sum(nn) ## Working 'totrs'
  ns <- max(strata)
  nstra <- c(0, cumsum(table(strata)))
  
  counts <- .C("sizes",
               as.integer(ns),
               as.integer(nn), 
               as.double(enter),
               as.double(exit),
               as.integer(event),
               ##
               antrs = integer(ns),
               as.integer(nstra),
               risktimes = double(w.totrs),
               ##
               n.events = integer(w.totrs),
               size = integer(w.totrs),
               totrs = integer(1),
               ## DUP = FALSE,
               PACKAGE = "eha"
               )

  counts$risktimes <- counts$risktimes[1:counts$totrs]
  counts$n.events <- counts$n.events[1:counts$totrs]
  counts$size <- counts$size[1:counts$totrs]

  totsize <- sum(counts$size)
  if (totsize >= 2^31 & members) stop("Too large  risk sets.") 
  totevents <- sum(counts$n.events)

  Eventset <- NULL
  Riskset <- NULL
  sample_fraction <- NULL
  
  if (members){
      res <- .C("risk_get",
                as.integer(max.survs),
                as.integer(nn),
                as.integer(ns),
                ##
                as.double(enter),
                as.double(exit),
                as.integer(event),
                ##
                as.integer(nstra),
                as.integer(length(nstra)),
                ##
                new.totrs = integer(1),  ## If sampling...
                ##
                as.integer(counts$antrs),
                as.integer(counts$n.events),
                size = as.integer(counts$size), ## If sampling...
                as.double(counts$risktimes),
                eventset = integer(totevents),
                riskset = integer(totsize),
                ## DUP = FALSE,
                PACKAGE = "eha"
                )
      Size <- res$size ## Previously out-commented; why??
      N <- counts$size - counts$n.events
      sample_fraction <- numeric(length(Size))
      pos <- N > 0
      sample_fraction[pos] <- (Size[pos] - counts$n.events[pos]) /
          (counts$size[pos] - counts$n.events[pos])
      sample_fraction[!pos] <- 1 # No survivors!
      Eventset <- ord[res$eventset]
      Riskset <- ord[res$riskset[1:res$new.totrs]]
  }

  rs <- list(ns = ns,
             antrs = counts$antrs,
             risktimes = counts$risktimes,
             n.events = counts$n.events,
             size = counts$size
             ##size = Size,
             ##eventset = Eventset,
             ##riskset = Riskset,
             ##sample_fraction = sample_fraction
             )
  if (members){
      rs$size <- Size
      rs$eventset <- Eventset
      rs$riskset <- Riskset
      rs$sample_fraction <- sample_fraction
      if (collate_sets){
          out <- vector(mode = "list", length = sum(rs$antrs))
          rstart <- 1
          estart <- 1
          start <- 1
          for (strata in 1:rs$ns){
        
              for (i in 1:rs$antrs[strata]){
            
                  ##cat("stratum = ", strata, "rs = ", i, "\n")
                  out[[start]] <-
                      list(stratum = strata,
                           risktime = rs$risktimes[rstart],
                           events =
                               rs$eventset[estart:(estart +
                                                   rs$n.events[start] - 1)],
                           risk = rs$riskset[rstart:(rstart +
                                                rs$size[start] - 1)],
                           sfrac = rs$sample_fraction[start])
                  estart <- estart + rs$n.events[start]
                  rstart <- rstart + rs$size[start]
                  start <- start + 1
              }
          }
          rs <- out
      }
  }
  class(rs) <- "risksets"
  
  invisible(rs)
}
