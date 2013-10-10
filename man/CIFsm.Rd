\name{CIFsm}
\alias{CIFsm}
\title{
Cumulative incidence function estimate and weighted summary statistics
}
\description{
Estimate the cumulative incidence function for cause of interest in two-sample study. Provide the weighted summary statistics based on given method and weight. 
}
\usage{
CIFsm(ds,method="dif",pp = 0,qq = 0,conf.bd=T,n.sim=500)
}

\arguments{
  \item{ds}{\code{ds} is a dataset contains the time, cause of event and group. 
  For cause, 0 means censoring, 1 is the cause of event, 2 is all other causes. Two groups need to be coded as 1 and 2.} 
  \item{method}{\code{method} can be chosen from "dif"=risk difference, "rr"=risk ratio and "or"=odds ratio}
  \item{pp}{first parameter of weight function}
  \item{qq}{second parameter of weight function}
  \item{conf.bd}{logical; if TRUE, create confidence band cut point. Set to FALSE if the confidence band is not needed, 
which reduces the computational time}
  \item{n.sim}{number of simulations used in creating confidence band; will be ineffective if \code{conf.bd} is FALSE}
}

\details{
The estimates and summary statistics are described in Zhang and Fine (2008).
}
\value{
\item{sample}{Total sample size from both groups}
\item{used}{Sample used in analysis. Subject with missing value in any of the three variables (time, cause or group) will be excluded from analysis}
\item{size}{Sample size for each group}
\item{njp}{Total number of unique event time points in two groups}
\item{tjp}{Unique event time points from both groups}
\item{ny1}{Number of subject at risk in group 1}
\item{f1}{Estimate of cumulative incidence function for group 1}
\item{f1.se}{Standard error of cumulative incidence function for group 1}
\item{ny2}{Number of subject at risk in group 2}
\item{f2}{Estimate of cumulative incidence function for group 2}
\item{f2.se}{Standard error of cumulative incidence function for group 2}
\item{dif}{Estimate of difference in cumulative incidence (risk difference) between the 2 groups}
\item{dif.se}{Standard error of the risk difference}
\item{dif.pv}{P-value of the risk difference}
\item{rr}{Risk ratio of the cumulative incidence between the 2 groups}
\item{rr.se}{Standard error of the risk ratio}
\item{rr.pv}{P-value of the risk ratio}
\item{or}{Odds ratio of the cumulative incidence between the 2 groups}
\item{or.se}{Standard error of the odds ratio}
\item{or.pv}{P-value of the odds ratio}
\item{cbcut}{95\% confidence band cut point of risk difference, risk ratio and odds ratio based on simulation method.}
\item{method}{Method used for the weighted summary statistics}
\item{weight}{Weight used for the weighted summary statistics}
\item{region}{Time range of the data}
\item{nbd}{Index of beginning and end of the data. For internal use only}
\item{ave}{Time integrated weighted summary statistics}
\item{avese}{Standard error of the time integrated weighted summary statistics}
\item{ci95}{95\% confidence interval of the time integrated weighted summary statistics}
\item{avepval}{P-value of the time integrated weighted summary statistics}
\item{wt}{Weight assigned at each unique event time point}
}
\references{
M.J. Zhang and J.P. Fine. Summarizing difference in cumulative incidence functions. 
Statistics in Medicine, 27:4939-4949, 2008.

J. Li, M.J. Zhang and J. Le-Rademacher. Weighted Summary of two Cumulative Incidence
Functions with R-CIFsmry Package. Computer Methods and Programs in Biomedicine[Submitted].
}

\author{
Jianing Li
}

\examples{
library(CIFsmry)
data(sim.dat)
out <- CIFsm(sim.dat,pp=0,qq=0)
out$avepval
plot(out$tjp,out$f1,type="s",ylim=c(0,1),yaxt="n",xaxt="n",xlab="Time",ylab="CIF",
     lty=1,lwd=1)
points(out$tjp,out$f2,type="s",lty=1,lwd=3)
axis(1,at=seq(0,6,by=1),cex=0.6)
axis(2,at=seq(0,1,by=0.2),cex=0.6)
legend("bottomright",c("1","2"),title="group",lty=1,lwd=c(1,3))

out10 <- CIFsm(sim.dat,pp=1,qq=0)
out10$avepval
}

\keyword{cumulative incidence function}
\keyword{summary statistics}
\keyword{weight function}

