\name{swepop}
\alias{swepop}
\docType{data}
\title{
Swedish population data, 1968-2019.
}
\description{
A data frame containing data on the population size by sex, age and
year, Sweden 1968-2019.
}
\usage{data(swepop)}
\format{
  A data frame with 5 variables and 10504 observations.
  \describe{
    \item{\code{age}}{Numerical with integer values 0-100, representing achieved
    age in years during the actual calendar year.The highest value, 100, 
    represents ages 100 and above.}
    \item{\code{sex}}{A factor with two levels, "women" and "men".}
    \item{\code{year}}{Calendar year.}
    \item{\code{pop}}{Average population by age, sex, and year.}
    \item{id}{Created by the \code{reshape} procedure, see Details.}
  }
}
\details{
Data are downloaded from Statistics Sweden in the form of a csv file and
converted to a data frame. Variable names are translated from Swedish,
and some of them are coverted to factors. The variable \code{pop} contains
the average population by sex and age, calculated by taking the mean
value of the population size at December 31 the previous year and
December 31 the current year. The original data contain the sizes at the
end of each year. The original data set is in wide format and converted to 
long format by \code{reshape}.
}
\source{
Statistics Sweden, \url{https://scb.se}.
}

\seealso{\code{\link{swedeaths}}}
  
\examples{
summary(swepop)
## maybe str(swepop) ...
}
\keyword{datasets}
