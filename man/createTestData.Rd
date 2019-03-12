\name{createTestData}
\alias{createTestData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a test data set for the example
}
\description{
Helper function for example that creates a test data set in the proper format
}
\usage{
createTestData(times = c(0:11), transition = 5, rx1 = 0.5, rx2 = 1, slope = 0.5,
error=0.3, n = c(50, 50, 50))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{Observation times
%%     ~~Describe \code{times} here~~
}
  \item{transition}{Time where phase I ends and Phase II begins
%%     ~~Describe \code{transition} here~~
}
  \item{rx1}{Drug effect on slope of outcome during phase 1.
%%     ~~Describe \code{rx1} here~~
}
  \item{rx2}{Drug effect on slope of outcome during phase 2.
%%     ~~Describe \code{rx2} here~~
}
  \item{slope}{Slope of outcome.
%%     ~~Describe \code{slope} here~~
}
  \item{error}{Standard Deviation of the error.}
  \item{n}{Sample sizes of the three groups, Placebo-Placebo, Placebo-Drug, Drug
%%     ~~Describe \code{n} here~~
}
}
\details{The model used is
that the outcome is equal to (\code{slope}) times \code{times} for patients in the Placebo cohort with the addition of \code{rx1} for the Drug cohort, during the first phase, and is equal to (\code{slope}+\code{rx2}) times (\code{times}-\code{transition}) for the Drug group in the second phase, and \code{slope} for the Placebo group in the second phase. A random noise is added to each observation with standard deviation \code{error}
}
\value{
A dataframe is returned with variables:
\item{ID}{Patient ID}
\item{rx}{Integer 1,2,3 indicating whether the patient is in the placebo-placebo, placebo-drug or drug groups}
\item{transition}{Time that phase 1 ends and phase 2 starts}
\item{time}{The observation time}
\item{y}{The value of the outcome variable}
 }
\references{
Fava, M., Evins, A. E., Dorer, D. J., and Schoenfeld, D. A. (2003). The problem of the placebo response in clinical trials for psychiatric disorders: culprits, possible remedies, and a novel study design approach. \emph{Psychotherapy and Psychosomatics}, \bold{72},3, 115--127.

Chen, Y. F., Yang, Y., Hung, H. J., & Wang, S. J. (2011). Evaluation of performance of some enrichment designs dealing with high placebo response in psychiatric clinical trials. \emph{Contemporary clinical trials}, \bold{32},4, 592-604.
}
\author{
David A. Schoenfeld
}

\examples{
createTestData(times=0:3,n=c(3,3,3),transition=5,rx1=0,rx2=0)
  }

