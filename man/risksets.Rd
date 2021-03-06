% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/risksets.R
\name{risksets}
\alias{risksets}
\title{Finds the compositions and sizes of risk sets}
\usage{
risksets(
  x,
  strata = NULL,
  max.survs = NULL,
  members = TRUE,
  collate_sets = FALSE
)
}
\arguments{
\item{x}{A \code{Surv} object.}

\item{strata}{Stratum indicator.}

\item{max.survs}{Maximum number of survivors in each risk set. If smaller
than the 'natural number', survivors are sampled from the present ones. No
sampling if missing.}

\item{members}{If TRUE, all members of all risk sets are listed in the
resulting list, see below.}

\item{collate_sets}{logical. If TRUE, group information by
risk sets in a list. Only if \code{members = TRUE}.}
}
\value{
A list with components (if \code{collate_sets = FALSE})
\item{antrs}{No. of risk sets in each
stratum. The number of strata is given by \code{length(antrs)}.}
\item{risktimes}{Ordered distinct failure time points.}
\item{eventset}{If
'members' is TRUE, a vector of pointers to events in each risk set, else
NULL.}
\item{riskset}{If 'members' is TRUE, a vector of pointers to the
members of the risk sets, in order. The 'n.events' first are the events. If
'members' is FALSE, 'riskset' is NULL.}
\item{size}{The sizes of the risk
sets.}
\item{n.events}{The number of events in each risk set.}
\item{sample_fraction}{If 'members' is TRUE, the sampling fraction of
survivors in each risk set.
}
}
\description{
Focus is on the risk set composition just prior to a failure.
}
\details{
If the input argument max.survs is left alone, all survivors are accounted
for in all risk sets.
}
\note{
Can be used to "sample the risk sets".
}
\examples{

 enter <- c(0, 1, 0, 0)
 exit <- c(1, 2, 3, 4)
 event <- c(1, 1, 1, 0)
 risksets(Surv(enter, exit, event))

}
\seealso{
\code{\link{table.events}}, \code{\link{coxreg}}.
}
\author{
Göran Broström
}
\keyword{survival}
