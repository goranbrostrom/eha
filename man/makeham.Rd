% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Makeham.R
\name{makeham}
\alias{makeham}
\alias{pmakeham}
\alias{dmakeham}
\alias{qmakeham}
\alias{hmakeham}
\alias{Hmakeham}
\alias{rmakeham}
\title{The Gompertz-Makeham Distribution}
\usage{
dmakeham(x, shape = c(1, 1), scale = 1, log = FALSE)
pmakeham(q, shape = c(1, 1), scale = 1, lower.tail = TRUE, log.p = FALSE)
qmakeham(p, shape = c(1, 1), scale = 1, lower.tail = TRUE, log.p = FALSE)
hmakeham(x, shape = c(1, 1), scale = 1, log = FALSE)
Hmakeham(x, shape = c(1, 1), scale = 1, log.p = FALSE)
rmakeham(n, shape = c(1, 1), scale = 1)
}
\arguments{
\item{shape}{A vector, default value c(1, 1).}

\item{scale}{defaulting to 1.}

\item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P(X \le
x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.}

\item{x, q}{vector of quantiles.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is
taken to be the number required.}

\item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\value{
\code{dmakeham} gives the density, \code{pmakeham} gives the distribution
function, \code{qmakeham} gives the quantile function, \code{hmakeham} gives the
hazard function, \code{Hmakeham} gives the cumulative hazard function, and
\code{rmakeham} generates random deviates.

Invalid arguments will result in return value \code{NaN}, with a warning.
}
\description{
Density, distribution function, quantile function, hazard function,
cumulative hazard function, and random generation for the Gompertz-Makeham 
distribution with parameters \code{shape} and \code{scale}.
}
\details{
The Gompertz-Makeham distribution with \code{shape} parameter \eqn{a} and \code{scale}
parameter \eqn{\sigma}{b} has hazard function given by 
\deqn{h(x) = a[2] + a[1] \exp(x/\sigma)}{h(x) = a[2] + a[1] exp(x/b)}
for \eqn{x \ge 0}{x >= 0}.
}
\keyword{distribution}
