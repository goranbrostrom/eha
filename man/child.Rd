\name{child}
\alias{child}
\docType{data}
\title{Child mortality, Skellefteå, Sweeden 1850--1900.}
\description{Children born in Skellefteå, Sweden, 1850-1884, are followed 
fifteen years or until death or out-migration.}

\usage{data(child)}
\format{
  A data frame with 26855 children born 1850-1884.
  \describe{
    \item{\code{id}}{An identification number.}
    \item{\code{m.id}}{Mother's id.}
    \item{\code{sex}}{Sex.}
    \item{\code{socBranch}}{Working branch of family (father).}
    \item{\code{birthdate}}{Birthdate.}
    \item{\code{enter}}{Start age of follow-up, always zero.}
    \item{\code{exit}}{Age of departure, either by death or emigration.}
    \item{\code{event}}{Type of departure, death = 1, right censoring = 0.}
    \item{\code{illeg}}{Born out of marriage ("illegitimate")?}
    \item{\code{m.age}}{Mother's age.}
  }
}
\details{The Skellefteå region is a large region in the northern part of Sweden.
}
\source{
Data originate from the Centre for Demographic and Ageing Research, 
Umeå University, Umeå, Sweden,
\url{https://www.umu.se/en/centre-for-demographic-and-ageing-research/}.}

\examples{
fit <- coxreg(Surv(enter, exit, event) ~ sex + socBranch, data = child, coxph = TRUE)
summary(fit)
}
\keyword{datasets}
