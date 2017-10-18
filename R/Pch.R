#' The Piecewise Constant Hazards distribution.
#'
#' Density, distribution function, quantile function, hazard function,
#' cumulative hazard function, and random generation for the Piecewice
#' Constant Hazards (pch) distribution.
#' 
#' @details 
#' The pch distribution has a hazard function that is piecewise constant
#' on intervals defined by cutpoints 
#' \deqn{0 < c_1 < \cdots < c_n < \infty, n \ge 0}{0 < c_1 < ... < c_n < \infty, n \ge 0}.
#' If \code{n = 0}, this reduces to an exponential distribution.
#' @name Pch
#' @aliases pch ppch dpch hpch Hpch qpch rpch
#' @usage 
#' ppch(q, cuts, levels, lower.tail = TRUE, log.p = FALSE)
#' dpch(x, cuts, levels, log = FALSE)
#' hpch(x, cuts, levels, log = FALSE)
#' Hpch(x, cuts, levels, log.p = FALSE)
#' qpch(p, cuts, levels, lower.tail = TRUE, log.p = FALSE)
#' rpch(n, cuts, levels)
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param cuts Vector of cut points defining the intervals where the hazard function
#' is constant.
#' @param levels Vector of levels (values of the hazard function). 
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken 
#' to be the number required.
#' @return \code{dpch} gives the density, 
#' \code{ppch} gives the distribution function, 
#' \code{qpch} gives the quantile function, 
#' \code{hpch} gives the hazard function, 
#' \code{Hpch} gives the cumulative hazard function, and
#' \code{rpch} generates random deviates.
#' @note the parameter \code{levels} must have length at least 2, and the 
#' number of cut points must be one less than the number of levels. q
#' @keywords distribution
#' @export
ppch <- function(q, cuts, levels, lower.tail = TRUE, log.p = FALSE){
    y <- Hpch(q, cuts, levels)
    if (log.p){
        if (lower.tail){
            y <- log(-expm1(-y)) ## log(-expm1(-y)) = log(1 - exp(-y))?!
        }else{
            y <- -y
        }
    }else{
        if (lower.tail){
            y <- -expm1(-y)
        }else{
            y <- exp(-y)
        }
    }
    y
}

#' @export        
dpch <- function(x, cuts, levels, log = FALSE){
    y <- hpch(x, cuts, levels) * ppch(x, cuts, levels, lower.tail = FALSE)
    if (log) y <- log(y)
    y
}

#' @export
hpch <- function(x, cuts, levels, log = FALSE){
    cuts <- sort(unique(cuts))
    p <- length(levels)
    if (length(cuts) != (p - 1))stop("Must be one more level than cut.")
    if (any(cuts <= 0)) stop("all cuts must be positive")
    if (any(x < 0)) stop("x must be all positive.")
    y <- numeric(length(x))
    cuts <- c(0, cuts, Inf)
    y[(cuts[1] <= x) & (x <= cuts[1 + 1])] <- levels[1]
    if (p > 1.5){
        for (i in 2:p){
            y[(cuts[i] < x) & (x <= cuts[i + 1])] <- levels[i]
        }
    }
    if (log) y <- log(y)
    y
}

#' @export    
Hpch <- function(x, cuts, levels, log.p = FALSE){
    cuts <- sort(unique(cuts))
    p <- length(levels)
    if (length(cuts) != (p - 1))stop("Must be one more level than cut.")
    if (any(cuts <= 0)) stop("all cuts must be positive")
    if (any(x < 0)) stop("x must be all positive.")
    y <- numeric(length(x))
    cuts <- c(0, cuts, Inf)
    who <- (cuts[1] <= x) & (x <= cuts[1 + 1])
    if (sum(who)){
        z <- x[who]
        y[who] <- levels[1] * z
    }
    su <- levels[1] * cuts[2]
    if (p > 1.5){
        for (i in 2:p){
            who <- (cuts[i] < x) & (x <= cuts[i + 1])
            if (sum(who)){
                y[who] <- su + levels[i] * (x[who] - cuts[i])
            }
            su <- su + levels[i] * (cuts[i + 1] - cuts[i])
        }
    }
    if (log.p) y <- log(y)
    y
}

#' @export            
qpch <- function(p, cuts, levels, lower.tail = TRUE, log.p = FALSE){
    if (log.p) p <- exp(p)
    if (any(p >= 1)) stop("p must be < 1") 
    if (any(p <= 0)) stop("p must be > 0") 
    if (!lower.tail) p <- 1 - p
    n <- length(p)
    y <- numeric(n)
    f <- function(q, x){
        ppch(q, cuts, levels) - x
    }
    for (i in 1:n){
        y[i] <- uniroot(f, interval = c(0, 2000), x = p[i])$root
    }
    y
}

#' @export
rpch <- function(n, cuts, levels){
    x <- runif(n)
    qpch(x, cuts, levels) # ', cuts, levels' added in 2.4-4.
}
