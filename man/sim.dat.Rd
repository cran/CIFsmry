\name{sim.dat}
\alias{sim.dat}
\non_function{}
\title{Simulated competing risk data}
\description{
  Simulated data with 100 subjects and 3 variates: time, cause and group.
}
\format{
The data has 100 rows and 3 columns.
\describe{
\item{time}{a numeric vector. Survival time.}
\item{cause}{a numeric vector code. 
Survival status. 1: failure from the cause of interest; 2: failure from other causes; 0: censored.}
\item{group}{a numeric vector code. 1: group 1; 2: group 2.}
}
}
\source{
Simulated data  
}
\examples{
data(sim.dat)
names(sim.dat)
table(sim.dat$cause)
table(sim.dat$group)
tapply(sim.dat$time,sim.dat$group,summary)
}
\keyword{datasets}
