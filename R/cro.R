#' Creates a minimal representation of a data frame.
#' 
#' Given a data frame with a defined response variable, this function creates a
#' unique representation of the covariates in the data frame, vector (matrix)
#' of responses, and a pointer vector, connecting the responses with the
#' corresponding covariates.
#' 
#' The rows in the data frame are converted to text strings with \code{paste}
#' and compared with \code{match}.
#' 
#' @param dat A data frame
#' @param response The column(s) where the response resides.
#' @return A list with components \item{y}{The response.} \item{covar}{A data
#' frame with unique rows of covariates.} \item{keys}{Pointers from \code{y} to
#' \code{covar}, connecting each response with its covariate vector.}
#' @note This function is based on suggestions by Anne York and Brian Ripley.
#' @author Göran Broström
#' @seealso \code{\link{match}}, \code{\link{paste}}
#' @keywords manip
#' @examples
#' 
#' dat <- data.frame(y = c(1.1, 2.3, 0.7), x1 = c(1, 0, 1), x2 = c(0, 1, 0))
#' cro(dat)
#' 
#' @export
cro <- function(dat, response = 1){
  covar <- unique(dat[, -response, drop = FALSE])
  dat.keys <-
    match(do.call("paste", c(dat[, -response, drop = FALSE], sep="\001")),
          do.call("paste", c(covar,  sep="\001")))

  list(y = dat[, response],
       covar = covar,
       keys = dat.keys)
}
