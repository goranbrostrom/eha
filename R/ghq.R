#' Gauss-Hermite
#' 
#' Calculates the zeros and weights needed for Gauss-Hermite quadrature.
#' 
#' Based on a Fortran 66 subroutine written by professor Jianming Jin.
#' 
#' @param n.points Number of points.
#' @param modified Multiply by exp(zeros**2)? Default is TRUE.
#' @return A list vith components \item{zeros}{The zeros (abscissas).}
#' \item{weights}{The weights}
#' @note The code is modified to suit the purpose of glmmML, with the
#' permission of professor Jin.
#' @author Jianming Jin, Univ. of Illinois, Urbana-Campaign
#' @seealso \code{\link{glmmML}}
#' @references Gauss-Hermite
#' @keywords math
#' @examples
#' 
#' ghq(15, FALSE)
#' 
#' @export ghq
ghq <- function(n.points = 1, modified = TRUE){
    weights <- numeric(n.points)
    zeros <- numeric(n.points)
    res <- .Fortran("ghq",
                    as.integer(n.points),
                    zeros = as.double(zeros),
                    weights = as.double(weights),
                    as.logical(modified),
                    PACKAGE = "eha")
    list(weights = res$weights, zeros = res$zeros)
}
