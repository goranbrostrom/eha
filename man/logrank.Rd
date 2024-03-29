% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logrank.R
\name{logrank}
\alias{logrank}
\title{The Log-rank test}
\usage{
logrank(Y, group, data = parent.frame())
}
\arguments{
\item{Y}{a survival object as returned by the \code{\link{Surv}} function.}

\item{group}{defines the groups to be compared. Coerced to a factor.}

\item{data}{a data.frame in which to interpret the variables.}
}
\value{
A list of class \code{logrank} with components
\item{test.statistic}{The logrank (score) test statistic.}
\item{df}{The degrees of freedom of the test statistic.}
\item{p.value}{The p value of the test.}
\item{hazards}{A list of two-column matrices, describing event times and 
corresponding hazard atoms in each stratum (class 'hazdata').}
\item{call}{The call}
}
\description{
Performs the log-rank test on survival data, possibly stratified.
}
\note{
The test is performed by fitting a Cox regression model and reporting
its \code{score test}. With tied data, this might be slightly different from
the true logrank test, but the difference is unimportant in practice.
}
\examples{
fit <- logrank(Y = Surv(enter, exit, event), group = civ, 
data = oldmort[oldmort$region == "town", ])
fit
}
\seealso{
\code{\link{coxreg}}, \code{\link{print.logrank}}.
}
\author{
Göran Broström
}
\keyword{survival}
