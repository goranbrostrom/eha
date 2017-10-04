#' Check the integrity of survival data.
#' 
#' Check that exit occurs after enter, that spells from an individual do not
#' overlap, and that each individual experiences at most one event.
#' 
#' Interval lengths must be strictly positive.
#' 
#' @param enter Left truncation time.
#' @param exit Time of exit.
#' @param event Indicator of event. Zero means 'no event'.
#' @param id Identification of individuals.
#' @param eps The smallest allowed spell length or overlap.
#' @return A vector of id's for the insane individuals. Of zero length if no
#' errors.
#' @author Göran Broström
#' @seealso \code{\link{join.spells}}, \code{\link{coxreg}},
#' \code{\link{aftreg}}
#' @keywords manip survival
#' @examples
#' 
#' xx <- data.frame(enter = c(0, 1), exit = c(1.5, 3), event = c(0, 1), id =
#' c(1,1))
#' check.surv(xx$enter, xx$exit, xx$event, xx$id)
#' 
#' @export
check.surv <- function(enter, exit, event, id = NULL, eps = 1.e-8){
    ## The '.Fortran' version.
    ##########################
    n <- length(enter)
    if (length(exit) != n)stop("Length mismatch (enter/exit)")
    if (length(event) != n)stop("Length mismatch (enter/event)")
    if(!is.null(id)) if (length(id) != n)stop("Length mismatch (enter/id)")

    ## If no id (or one record per id).
    if (is.null(id) || (length(unique(id)) == n)) return(all(enter < exit))

    ## Now, id is set; let's sort data:
    #id <- factor(id)
    n.ind <- length(unique(id))
    ord <- order(id, enter)
    id <- id[ord]
    enter <- enter[ord]
    exit <- exit[ord]
    event <- as.logical(event[ord])

    id <- factor(id) 
    id.size <- table(id)

    xx <- .Fortran("chek",
                   as.integer(n),
                   as.integer(n.ind),
                   as.integer(id.size),    ## length = n.ind
                   as.double(enter),       ## length = n
                   as.double(exit),        ## length = n
                   as.integer(event),      ## length = n
                   as.double(eps),
                   sane = integer(n.ind)   ## boolean; TRUE: good individual
                   )

    bad.id <- levels(id)[xx$sane == 0]
    bad.id
}
