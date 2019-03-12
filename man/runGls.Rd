\name{runGls}
\alias{runGls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Runs a constrained longitudinal data analysis model for use with \code{SPCDcontinuous}
}
\description{ The function \code{SPCDFcontinuous} requires an analysis function which estimates the treatment effect and it's standard error for each phase of the study. This function uses the so called cLDA model to analyze the data.  In the situation where their are only pre-test and post-test measurements this is roughly equivalent to an analysis of covariance.  This function can either be used as is or be used as a model for a user supplied function.
}
\usage{
runGls(combData, times = "time")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{combData}{This is the input file data file, it needs to have fields \code{ID} which is the patient identifier,
\code{time} which is the observation time, \code{y} the outcome variable, \code{trt} a treatment code which is either
\code{Placebo} or \code{Active}
}
  \item{times}{
This is the time variable. The set of times for each patient needs to be the same(except for some patients times may be missing).  If they are different then the \code{runSlopes} function might work.
}
}

\value{
A numerical vector of length 3 which gives the treatment effect, it's standard are and the z-score which is the ratio of the treatment effect to it's standard error. The function \code{SPCDcontinuous} will use the number of patients as the degrees of freedom to compute p-values.
}
\references{
Liang, K. Y. and Zeger, S. (2000). Longitudinal data analysis of continuous and discrete responses for pre-post designs. \emph{Sankhya: The Indian Journal of Statistics, Series B} \bold{62}, 134–148.
}
\author{
David A. Schoenfeld dschoenfeld@mgh.harvard.edu
}

\seealso{
\code{\link{SPCDcontinuous}}
}
\examples{

data1=createTestData(times=c(0,5,11),n=c(50,50,50),transition=5,rx1=0,rx2=0)
data1=placeboNonResponder(data1,k=0)
#Example of use of function outside of SPCDcontinous
data1$trt=ifelse(data1$rx==3,'Active','Placebo')
runGls(data1,time='time')
}
\keyword{htest}% use one of  RShowDoc("KEYWORDS")

