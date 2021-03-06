% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phreg.fit.R
\name{phreg.fit}
\alias{phreg.fit}
\title{Parametric proportional hazards regression}
\usage{
phreg.fit(X, Y, dist, strata, offset, init, shape, control)
}
\arguments{
\item{X}{The design (covariate) matrix.}

\item{Y}{A survival object, the response.}

\item{dist}{Which baseline distribution?}

\item{strata}{A stratum variable.}

\item{offset}{Offset.}

\item{init}{Initial regression parameter values.}

\item{shape}{If positive, a fixed value of the shape parameter in the
distribution. Otherwise, the shape is estimated.}

\item{control}{Controls convergence and output.}
}
\value{
\item{coefficients}{Estimated regression coefficients plus estimated
scale and shape coefficients, sorted by strata, if present.}
\item{var}{Variance-covariance matrix} \item{loglik}{Vector of length 2. The
first component is the maximized loglihood with only scale and shape in the
model, the second the final maximum.} \item{score}{Score test statistic at
initial values} \item{linear.predictors}{Linear predictors for each
interval.} \item{means}{Means of the covariates} \item{conver}{TRUE if
convergence} \item{fail}{TRUE if failure} \item{iter}{Number of
Newton-Raphson iterates.} \item{n.strata}{The number of strata in the data.}
}
\description{
This function is called by \code{\link{phreg}}, but it can also be directly
called by a user.
}
\details{
See \code{\link{phreg}} for more detail.
}
\seealso{
\code{\link{phreg}}
}
\author{
Göran Broström
}
\keyword{regression}
\keyword{survival}
