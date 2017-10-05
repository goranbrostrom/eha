\name{oldmort}
\alias{oldmort}
\docType{data}
\title{
Old age mortality, Sundsvall, Sweden, 1860-1880.
}
\description{
The data consists of old age life histories from 1 January 1860 to 31 december
1880, 21 years. Only (parts of) life histories above age 60 is
considered. 
}
\usage{data(oldmort)}
\format{
  A data frame with 6508 observations from 4603 persons on the following 13 variables.
  \describe{
    \item{\code{id}}{Identification number.}
    \item{\code{enter}}{Start age for the interval.}
    \item{\code{exit}}{Stop age for the interval.}
    \item{\code{event}}{Indicator of death; equals \code{TRUE} if the person died
      at the end of the interval, \code{FALSE} otherwise.}
    \item{\code{birthdate}}{Birthdate as a real number (i.e.,
      "1765-06-27" is 1765.490).}
    \item{\code{m.id}}{Mother's identification number.}
    \item{\code{f.id}}{Father's identification number.}
    \item{\code{sex}}{Gender, a factor with levels \code{male} \code{female}}
    \item{\code{civ}}{Civil status, a factor with levels \code{unmarried} \code{married} \code{widow}}
    \item{\code{ses.50}}{Socio-economic status at age 50, a factor with levels \code{middle} \code{unknown} \code{upper} \code{farmer} \code{lower}}
    \item{\code{birthplace}}{a factor with levels \code{parish} \code{region} \code{remote}}
    \item{\code{imr.birth}}{Infant mortality rate at birth in the region
    of birth}
    \item{\code{region}}{Subregion of Sundsvall, a factor with levels \code{town} \code{industry} \code{rural}}
  }
}
\details{
The Sundsvall area in mid-Sweden was during the 19th century a fast
growing forest industry. At the end of the century, it was one of the
largest sawmill area in Europe. The town Sundsvall is fast growing part
of the region and center for the commerse.
}
\source{
The Demographic Data Base, Umeå University, Sweden.
}
\references{
Edvinsson, S. (2000). The Demographic Data Base at Umeå University: A
resource for historical studies. In Hall, McKaa, and Thorvaldsen (eds),
"Handbook of International Historical Microdata for Population
Research", Minnesota Population Center, Minneapolis.
}
\examples{
data(oldmort)
summary(oldmort)
## maybe str(oldmort) ; plot(oldmort) ...
}
\keyword{datasets}
