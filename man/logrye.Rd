\name{logrye}
\alias{logrye}
\docType{data}
\title{
Rye prices, Scania, southern Sweden, 1801-1894.
}
\description{
The data consists of yearly rye prices from 1801 to
1894. Logged and detrended, so the time series is supposed to measure
short term fluctuations in rye prices.
}
\usage{data(scania)}
\format{
  A data frame with 94 observations in two columns on the following 2 variables.
  \describe{
    \item{\code{year}}{The year the price is recorded.}
    \item{\code{foodprices}}{Detrended log rye prices.}
  }
}
\details{
The Scanian area in southern Sweden was during the 19th century a mainly
rural area. 
}
\source{
The Scanian Economic Demographic Database.
}
\references{Jörberg, L. (1972). A History of Prices in Sweden 1732-1914,
  CWK Gleerup, Lund.
}
\examples{
data(logrye)
summary(logrye)
}
\keyword{datasets}
