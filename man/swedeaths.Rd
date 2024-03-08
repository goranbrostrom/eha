\name{swedeaths}
\alias{swedeaths}
\docType{data}
\title{
Swedish death data, 1969-2020.
}
\description{
A data frame containing data on the number of deaths by sex, age and
year, Sweden 1969-2020.}
\usage{data(swedeaths)}
\format{
  A data frame with 5 variables and 10504 observations.
  \describe{
    \item{\code{age}}{Numerical with integer values 0-100, representing achieved
    age in years during the actual calendar year. The highest value, 100, 
    represents ages 100 and above.}
    \item{\code{sex}}{A factor with two levels, "women" and "men".}
    \item{\code{year}}{Calendar year.}
    \item{\code{deaths}}{Number of deaths by age, sex, and year.}
    \item{\code{id}}{ Created by the \code{reshape} procedure, see Details.}
  }
}
\details{
Data are downloaded from Statistics Sweden in the form of a csv file and
and in that process converted to a data frame. Variable names are
translated from Swedish, 
and some of them are coverted to factors. Each numeric column contains
the number of deaths by sex and age. The original data set is in wide
form and then converted to long format.
}
\source{
Statistics Sweden, \url{https://www.scb.se}.
}

\seealso{\code{\link{swepop}}, \code{\link{tpchreg}}}

\examples{
summary(swedeaths)
## maybe str(swedeaths) ...
}
\keyword{datasets}
