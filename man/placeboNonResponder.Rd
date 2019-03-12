\name{placeboNonResponder}
\alias{placeboNonResponder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Example of helper function for determining placebo nonresponse}
%%  ~~function to do ... ~~

\description{The input to \code{SPCDcontinuous} requires a logical field which says that the patient was not a placebo responder in the the first phase of the study. Usually this is determined by whether their outcome was greater than a specified value but could have other criteria.  This an example program which specifies non-responders as patients whose first phase results were greater than a constant.

}
\usage{
placeboNonResponder(data, ID = "ID", k = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data file to use to determine non-response.
}
  \item{ID}{
Text field indicating the name of the Patient ID variable.
}
  \item{k}{
A numeric constant.  A patient is a non-responder if their change from baseline was greater than k.
}
}
\value{A new data frame with a variable named \code{nonResponder}.
}
\references{
Chen, Y. F., Yang, Y., Hung, H. J., & Wang, S. J. (2011). Evaluation of performance of some enrichment designs dealing with high placebo response in psychiatric clinical trials. \emph{Contemporary clinical trials}, \bold{32},4, 592-604.
}
\author{
David A. Schoenfeld, dschoenfeld@mgh.harvard.edu
}
\note{
The choice of criteria should be specified in the protocol. One suggesting is to consider a non responder to be a patient who would still be eligible for the study.
}
\examples{
data1=createTestData(times=0:11,n=c(50,50,50),transition=5,rx1=0,rx2=0)
data1=placeboNonResponder(data1,k=0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
