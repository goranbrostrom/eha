\name{infants}
\alias{infants}
\docType{data}
\title{Infant mortality and maternal death, Sweeden 1821--1894.}
\description{Matched data on infant mortality, from seven parishes in
  Sweden, 1821--1894.}  

\usage{data(infants)}
\format{
  A data frame with 80 rows and five variables.
  \describe{
    \item{\code{stratum}}{Triplet No. Each triplet consist of one infant
      whose mother died (a case), and two controls, i.e, infants whose
      mother did 
      not die. Matched on covariates below.} 
    \item{\code{enter}}{Age (in days) of case when its mother died.}
    \item{\code{exit}}{Age (in days) at death or right censoring (at age
      365 days).}
    \item{\code{event}}{Follow-up ends with death (1) or right censoring
      (0).}
    \item{\code{mother}}{\code{dead} for cases, \code{alive} for controls.}
    \item{\code{age}}{Mother's age at infant's birth.}
    \item{\code{sex}}{The infant's sex.}
    \item{\code{parish}}{Birth parish, either Nedertornea or not
	Nedertornea.}
     \item{\code{civst}}{Civil status of mother, \code{married} or
       \code{unmarried}.}
   \item{\code{ses}}{Socio-economic status of mother, either farmer or
     not farmer.}
   \item{\code{year}}{Year of birth of the infant.}
  }
}
\details{
From 5641 first-born in seven Swedish parishes 1820-1895, from Fleninge in the
very south to Nedertorneå in the very north, those whose mother died
during their first year of life were selected, in all 35 infants. To
each of them, two controls were selected by matching on the given covariates.
} 
\source{
Data originate from The Demographic Data Base, Umeå University, Umeå, Sweden,
\url{http://www.cedar.umu.se}.}
\references{
Broström, G. (1987). The influence of mother's death on infant
mortality: A case study in matched data survival analysis. Scandinavian
Journal of Statistics 14, 113-123.}

\examples{
data(infants)
fit <- coxreg(Surv(enter, exit, event) ~ strata(stratum) + mother, data
= infants)
fit
fit.w <- phreg(Surv(enter, exit, event) ~ mother + parish + ses, data =
infants)
summary(fit.w) ## Weibull proportional hazards model.
}
\keyword{datasets}
