\name{scania}
\alias{scania}
\docType{data}
\title{
Old age mortality, Scania, southern Sweden, 1813-1894.
}
\description{
The data consists of old age life histories from 1 January 1813 to 31 december
1894. Only (parts of) life histories above age 50 is
considered. 
}
\usage{data(scania)}
\format{
  A data frame with 1931 observations from 1931 persons on the following 9 variables.
  \describe{
    \item{\code{id}}{Identification number (enumeration).}
    \item{\code{enter}}{Start age for the interval.}
    \item{\code{exit}}{Stop age for the interval.}
    \item{\code{event}}{Indicator of death;
      equals \code{TRUE} if the person died
      at the end of the interval, \code{FALSE} otherwise.}
    \item{\code{birthdate}}{Birthdate as a real number (i.e.,
      "1765-06-27" is 1765.490).}
    \item{\code{sex}}{Gender, a factor with levels \code{male}
      \code{female}.}
    \item{\code{parish}}{One of five parishes in Scania, coded 1, 2, 3,
	4, 5. Factor.}
    \item{\code{ses}}{Socio-economic status at age 50, a factor with
      levels \code{upper} and \code{lower}.}
    \item{\code{immigrant}}{a factor with levels \code{no} \code{region}
      and \code{yes}.}
  }
}
\details{
The Scanian area in southern Sweden was during the 19th century a mainly
rural area. 
}
\source{
The Scanian Economic Demographic Database, Lund University, Sweden.
}
\references{\url{http://www.ed.lu.se/databases}
}
\examples{
data(scania)
summary(scania)
}
\keyword{datasets}
