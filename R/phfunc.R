#' Loglihood function of a proportional hazards regression
#' 
#' Calculates minus the log likelihood function and its first and second order
#' derivatives for data from a Weibull regression model.
#' 
#' Note that the function returns log likelihood, score vector and minus
#' hessian, i.e. the observed information. The model is 
#' \deqn{S(t; p, \lambda, \beta, z) = S_0((t / \lambda)^p)^{e^(z \beta)}}{S(t; p, lambda, beta, z) = S_0((t / lambda)^p)^exp(z beta)}
#' 
#' @param beta Regression parameters
#' @param lambda The scale paramater
#' @param p The shape parameter
#' @param X The design (covariate) matrix.
#' @param Y The response, a survival object.
#' @param offset Offset.
#' @param ord ord = 0 means only loglihood, 1 means score vector as well, 2
#' loglihood, score and hessian.
#' @param pfixed Logical, if TRUE the shape parameter is regarded as a known
#' constant in the calculations, meaning that it is not cosidered in the
#' partial derivatives.
#' @param dist Which distribtion? The default is "weibull", with the
#' alternatives "loglogistic" and "lognormal".
#' @return A list with components \item{f}{The log likelihood. Present if
#' \code{ord >= 0}} \item{fp}{The score vector. Present if \code{ord >= 1}}
#' \item{fpp}{The negative of the hessian. Present if \code{ord >= 2}}
#' @author Göran Broström
#' @seealso \code{\link{phreg}}
#' @keywords survival distribution
#' @export phfunc
phfunc <- function(beta = NULL, lambda, p, X = NULL, Y,
                  offset = rep(0, length(Y)),
                  ord = 2, pfixed = FALSE, dist = "weibull"){

## Returns loglik, score, and information (=-fpp)
## For one stratum (only)!!

  if (ord < 0) return(NULL)

  if (dist == "weibull"){
    dis <- 0
  }else if(dist == "loglogistic"){
    dis <- 1
  }else if (dist == "lognormal"){
    dis <- 2
  }


    nn <- NROW(Y)
    if (NCOL(Y) == 2) Y <- cbind(rep(0, nn), Y)

    if (is.null(X)){
        stop("No X not implemented yet")
        if (FALSE){    #### NOT yet!!
            if (pfixed){
                bdim <- 1
                b <- log(lambda)
            }else{
                bdim <- 2
                b <- c(log(lambda), log(p))
            }

    ##        fit <- .Fortran("phfuncnull",
      ##                      as.integer(ord),
      ##                      as.integer(pfixed),
      ##                      as.double(p),
      ##                      as.integer(bdim),
      ##                      as.double(b),
      ##                      as.integer(nn),
      ##                      as.double(Y[, 1]),
      ##                      as.double(Y[, 2]),
      ##                      as.integer(Y[, 3]),
                                        #
      ##                      f = double(1),
      ##                      fp = double(bdim),
      ##                      fpp = double(bdim * bdim),
      ##                      ok = integer(1),
                            ## DUP = FALSE,
      ##                      PACKAGE = "eha"
      ##                      )
        }

    }else{
        mb <- NCOL(X)
        if (length(beta) != mb) stop("beta mis-specified!")
        if (pfixed){
          stop("p fixed not implemented yet")
            bdim <- mb + 1
            b <- c(beta, log(lambda))
        }else{
            bdim <- mb + 2
            b <- c(beta, log(lambda), log(p))
        }
        b <- beta
        alpha <- log(lambda)
        gamma <- log(p)

        d0 <- .C("loglik_ph",
                 as.integer(dis),
                 as.integer(mb),
                 as.double(b),
                 as.double(alpha),
                 as.double(gamma),
                 as.integer(nn),
                 as.double(t(X)),
                 as.double(Y[, 1]),
                 as.double(Y[, 2]),
                 as.integer(Y[, 3]),
                 as.double(offset),
                 f = double(1), ## Return value
                 ## DUP = FALSE,
                 PACKAGE = "eha")
        ret <- list(f = - d0$f)
        if (ord >= 1){
          d1 <- .C("d_loglik_ph",
                   as.integer(dis),
                   as.integer(mb),
                   as.double(b),
                   as.double(alpha),
                   as.double(gamma),
                   as.integer(nn),
                   as.double(t(X)),
                   as.double(Y[, 1]),
                   as.double(Y[, 2]),
                   as.integer(Y[, 3]),
                   as.double(offset),
                   fp = double(bdim), ## Return value
                   ## DUP = FALSE,
                   PACKAGE = "eha")
          ret$fp <- d1$fp
          if (ord >= 2){
            d2 <- .C("d2_loglik_ph",
                     as.integer(dis),
                     as.integer(mb),
                     as.double(b),
                     as.double(alpha),
                     as.double(gamma),
                     as.integer(nn),
                     as.double(t(X)),
                     as.double(Y[, 1]),
                     as.double(Y[, 2]),
                     as.integer(Y[, 3]),
                     as.double(offset),
                     fpp = double(bdim * bdim), ## Return value
                     ## DUP = FALSE,
                     PACKAGE = "eha")
            ret$fpp <- matrix(d2$fpp, ncol = bdim)
          }
        }
    }

  return(ret)
}

