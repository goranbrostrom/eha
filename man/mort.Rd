\name{mort}
\alias{mort}
\docType{data}
\title{Male mortality in ages 40-60, nineteenth century}
\description{
  Males born in the years 1800-1820 and surving at least 40 years  in
  the parish Skellefteå in northern 
  Sweden are followed from their fortieth birthday until death or the
  sixtieth birthday, whichever comes first.
}
\usage{data(mort)}
\format{
  A data frame with 2058 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{Personal identification number.} 
    \item{\code{enter}}{Start of duration. Measured in years since the
      fortieth birthday.}
    \item{\code{exit}}{End of duration. Measured in years since the
      fortieth birthday.}
    \item{\code{event}}{a logical vector indicating death at end of interval.}
    \item{\code{birthdate}}{The birthdate in decimal form.}
    \item{\code{ses}}{Socio-economic status, a factor with levels
      \code{lower}, \code{upper}} 
  }
}
\note{This data set is also known, and accessible,  as \code{male.mortality}}
\details{
  The interesting explanatory
  covariate is \code{ses} (socioeconomic status), which is a
  time-varying covariate. This explains why several individuals are
  representated by more than one record each. Left trucation and right
  censoring are introduced this way.
} 
\source{
  Data is coming from The Demographic Data Base, Umea University, Umeå, Sweden.
}
\references{
\url{http://www.cedar.umu.se}  
}
\examples{
data(mort)
fit <- coxreg(Surv(enter, exit, event) ~ ses, data = mort)
summary(fit)
}
\keyword{datasets}
