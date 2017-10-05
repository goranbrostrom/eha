\name{swe07}
\alias{swe07}
\docType{data}
\title{Swedish population and deaths in ages 61--80, 2007.}
\description{The Swedish population and No. of deaths by age and sex in
  the ages 61--80. Data from the year 2007. 
}
\usage{data(swe07)}
\format{
  A data frame with 80 rows and five variables.
  \describe{
    \item{\code{pop}}{Average population size during the year 2007 by
      age and sex.} 
    \item{\code{deaths}}{Number of deaths by age and sex during the year
    2007.}
    \item{\code{sex}}{Sex.}
    \item{\code{age}}{Age.}
    \item{\code{log.pop}}{The logarithm of the first variable,
      \code{pop}. Included for convenience, may be used as an offset in
      a Poisson regression.}
  }
}
\details{
The average population is calculated as the mean of the population 1
January 2007 and 1 January 2008.
} 
\source{
Data is taken from Statistics Sweden.
}
\references{
\url{http://www.scb.se}  
}
\examples{
data(swe07)
fit <- glm(deaths ~ offset(log.pop) + sex * as.factor(age), family = poisson, data = swe07)
drop1(fit, test = "Chisq")  ## Proportional hazards?
}
\keyword{datasets}
