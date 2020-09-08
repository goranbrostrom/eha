\name{fert}
\alias{fert}
\docType{data}
\title{Marital fertility nineteenth century}
\description{
Birth intervals for married women with at least one birth, 19th northern
Sweden 
}
\usage{data(fert)}
\format{
  A data frame with 12169 observations the lengths (in years) of birth
  intervals for 1859 married women with at least one birth. The first
  interval (\code{parity = 0}) is the interval from marriage to first birth.
  \describe{
    \item{\code{id}}{Personal identification number for mother.} 
    \item{\code{parity}}{Time order of birth interval for the present
      mother. The interval with \code{parity = 0} is the first, from
    marriage to first birth.}
    \item{\code{age}}{The age of mother at start of interval.}
    \item{\code{year}}{The calendar year at start of interval.}
    \item{\code{next.ivl}}{The length of the coming time interval.}
    \item{\code{event}}{An indicator for whether the \code{next.ivl}
      ends in a new birth (\code{event = 1}) or is right censored
      (\code{event = 0}). Censoring occurs when the woman ends her
      fertility period within her first marriage (marriage dissolution
      or reaching the age of 48).}
    \item{\code{prev.ivl}}{The length of the previous time interval. May
      be used as explanatory variable in a Cox regression of birth intervals.}
    \item{\code{ses}}{Socio-economic status, a factor with levels
      \code{lower}, \code{upper}, \code{farmer}, and \code{unknown}.}
    \item{\code{parish}}{The Skelleftea region consists of three
      parishes, Jorn, Norsjo, and Skelleftea.}
  }
}
\details{
  The data set contain clusters of dependent observations defined by
  mother's id. 
} 
\source{
  Data is coming from The Demographic Data Base, Umea University, Umea, Sweden.
}
\references{
\url{http://www.cedar.umu.se}  
}
\examples{
data(fert)
fit <- coxreg(Surv(next.ivl, event) ~ ses + prev.ivl, data = fert, subset =
(parity == 1))
summary(fit)
}
\keyword{datasets}
